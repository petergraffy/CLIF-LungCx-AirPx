# ================================================================================================
# ICU Post–Lung Cancer Resection Cohort Builder (CLIF) | PI: Peter Graffy
#
# Cohort: Adults (>=18) with lung cancer diagnosis present-on-admission (POA=1),
#         a qualifying lung resection procedure, and an ICU admission AFTER the procedure.
#
# Key linkage rules:
#   - Lung cancer dx: hospital_diagnosis where poa_present == 1 AND ICD code matches lung cancer codes
#   - Resection: patient_procedures where procedure_code matches resection codes
#   - ICU: adt where location_category == "ICU"
#   - Temporal: procedure_billed_dttm < ICU in_dttm
#   - Choose the closest ICU admission after the resection (per hospitalization)
#
# Outputs (output/run_[SITE]_[DATE]/):
#   cohort_lung_resection_icu.csv    : 1 row / hospitalization
#   exclusion_lung_resection_icu.csv : exclusions + first failing reason
#   flow_lung_resection_icu.csv      : flow counts
#
# ================================================================================================

suppressPackageStartupMessages({
  library(fst)
  library(here)
  library(tidyverse)
  library(arrow)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(data.table)
  library(readr)
  library(glue)
})

# ---------- Project config ----------
source("utils/config.R")
stopifnot(exists("config"))

repo        <- config$repo
site_name   <- config$site_name
tables_path <- config$tables_path
file_type   <- config$file_type

stopifnot(!is.null(repo), nzchar(repo))
stopifnot(!is.null(tables_path), nzchar(tables_path))

cat("Site Name:", site_name, "\n")
cat("Tables Path:", tables_path, "\n")
cat("File Type:", file_type, "\n")

tables_path <- normalizePath(tables_path, mustWork = TRUE)

# ---------- Cohort parameters ----------
START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")

ADULT_AGE_YEARS <- 18

# Require ICU admission occurs within this many hours AFTER the procedure (set to Inf to disable)
PROC_TO_ICU_MAX_H <- 48

# Index time choice for downstream recovery trajectories:
#   - "icu_in"  : t0 = ICU in time (recommended for post-op ICU recovery)
#   - "proc"    : t0 = procedure time
T0_ANCHOR <- "icu_in"

# Paths
CODES_PATH <- file.path(repo, "outlier-thresholds", "table.tsv")  # or wherever you put it
stopifnot(file.exists(CODES_PATH))

# ---------- Helpers ----------
safe_posix <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  suppressWarnings(as.POSIXct(x, tz = "UTC"))
}
safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "SITE" else x
}

# Normalize code strings (drop punctuation, uppercase)
norm_code <- function(x) {
  x <- toupper(as.character(x))
  x <- str_replace_all(x, "[^A-Z0-9]", "")
  x
}

# ICD often needs prefix matching (e.g., C34, C340, C341...)
code_matches_any_prefix <- function(code_vec, prefixes) {
  code_vec <- norm_code(code_vec)
  prefixes <- norm_code(prefixes)
  # true if any prefix matches start of code
  vapply(code_vec, function(cd) any(str_starts(cd, prefixes)), logical(1))
}

# ---------- File discovery / loading (same pattern as your REFER script) ----------
exts <- strsplit(file_type, "[/|,; ]+")[[1]]
exts <- exts[nzchar(exts)]
if (length(exts) == 0) exts <- c("csv","parquet","fst")
ext_pat <- paste0("\\.(", paste(unique(exts), collapse = "|"), ")$")

all_files <- list.files(
  path = tables_path,
  pattern = ext_pat,
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0) {
  stop("No files with extensions {", paste(exts, collapse = ", "), "} were found under: ", tables_path)
}

bn <- basename(all_files)
looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)
base_no_ext <- tools::file_path_sans_ext(tolower(bn))
base_norm <- ifelse(looks_clif, base_no_ext, paste0("clif_", base_no_ext))
found_map <- stats::setNames(all_files, base_norm)

required_raw <- c("patient","hospitalization","adt","hospital_diagnosis","patient_procedures","respiratory_support")
required_files <- paste0("clif_", required_raw)
missing <- setdiff(required_files, names(found_map))
if (length(missing) > 0) {
  cat("Detected CLIF-like files:\n"); print(sort(unique(names(found_map))))
  stop("Missing required tables: ", paste(missing, collapse = ", "))
}

clif_paths <- found_map[required_files]

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

clif_tables <- lapply(clif_paths, read_any)
names(clif_tables) <- required_files
cat("Loaded tables: ", paste(names(clif_tables), collapse = ", "), "\n")

get_min <- function(tbl_name, cols) {
  nm <- paste0("clif_", tbl_name)
  stopifnot(!is.null(clif_tables[[nm]]))
  out <- clif_tables[[nm]] %>% rename_with(tolower)
  cols_keep <- intersect(tolower(cols), names(out))
  out %>% dplyr::select(any_of(cols_keep))
}

# ---------- Load code table & build code vectors ----------
codes <- readr::read_tsv(CODES_PATH, show_col_types = FALSE) %>%
  rename_with(tolower)

# Expected columns (best-effort):
#   - code_system (e.g., "ICD10CM","ICD9CM","CPT","ICD10PCS","ICD9PCS")
#   - code
# If your column names differ, adjust here:
if (!all(c("code type","code") %in% names(codes))) {
  stop("Codes table must contain columns 'code_system' and 'code'. Found: ",
       paste(names(codes), collapse = ", "))
}

codes <- codes %>%
  mutate(
    code_system = toupper(str_replace_all(`code type`, "[^A-Z0-9]", "")),
    code = norm_code(code)
  )

# Split codes for dx vs procedure
lung_dx_prefixes <- codes %>%
  filter(code_system %in% c("ICD10CM","ICD9CM")) %>%
  pull(code) %>%
  unique()

lung_proc_codes <- codes %>%
  filter(code_system %in% c("CPT","ICD10PCS","ICD9PCS")) %>%
  pull(code) %>%
  unique()

cat("Loaded code counts:\n")
cat("  Lung dx prefixes (ICD9/10-CM): ", length(lung_dx_prefixes), "\n", sep="")
cat("  Lung resection procedure codes (CPT/PCS): ", length(lung_proc_codes), "\n", sep="")

# ---------- Minimal tables ----------
patient <- get_min("patient",
                   c("patient_id","birth_date","sex_category","race_category","ethnicity_category")) %>%
  mutate(birth_date = safe_date(birth_date))

hospitalization <- get_min("hospitalization",
                           c("patient_id","hospitalization_id","admission_dttm","discharge_dttm","age_at_admission",
                             "zipcode_nine_digit","zipcode_five_digit","census_tract","county_code")) %>%
  mutate(
    admission_dttm = safe_posix(admission_dttm),
    discharge_dttm = safe_posix(discharge_dttm)
  )

adt <- get_min("adt",
               c("hospitalization_id","in_dttm","out_dttm","location_category","location_type")) %>%
  mutate(
    in_dttm  = safe_posix(in_dttm),
    out_dttm = safe_posix(out_dttm)
  )

hospital_dx <- get_min("hospital_diagnosis",
                       c("hospitalization_id","diagnosis_code","diagnosis_code_format","poa_present","diagnosis_primary")) %>%
  mutate(
    icd_code = norm_code(diagnosis_code),
    poa_present = suppressWarnings(as.integer(poa_present))
  )

patient_proc <- get_min("patient_procedures",
                        c("hospitalization_id","procedure_code","procedure_code_format",
                          "procedure_billed_dttm")) %>%
  mutate(
    procedure_code = norm_code(procedure_code),
    procedure_code_type = toupper(str_replace_all(as.character(procedure_code_format), "[^A-Z0-9]", "")),
    procedure_billed_dttm = safe_posix(procedure_billed_dttm)
  ) %>%
  mutate(
    proc_time = coalesce(procedure_billed_dttm)
  )

# ---------- ICU bounds (first ICU in) ----------
icu_segments <- adt %>%
  mutate(is_icu = str_detect(tolower(coalesce(location_category, "")), "icu")) %>%
  filter(is_icu)

icu_bounds <- icu_segments %>%
  group_by(hospitalization_id) %>%
  summarize(
    first_icu_in = suppressWarnings(min(in_dttm, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
    last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours"))
  )

# ---------- Base ICU candidates (date range + demos + geo) ----------
base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  filter(!is.na(first_icu_in),
         first_icu_in >= START_DATE,
         first_icu_in <= END_DATE) %>%
  left_join(patient %>% select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
            by = "patient_id") %>%
  mutate(
    age_years = coalesce(
      suppressWarnings(as.numeric(age_at_admission)),
      ifelse(!is.na(birth_date),
             as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
    )
  ) %>%
  mutate(
    has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
    adult    = !is.na(age_years) & age_years >= ADULT_AGE_YEARS,
    has_geo  = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code)
  )

# ---------- Lung cancer POA identification ----------
lung_poa <- hospital_dx %>%
  filter(poa_present == 1) %>%
  filter(code_matches_any_prefix(icd_code, lung_dx_prefixes)) %>%
  distinct(hospitalization_id) %>%
  mutate(has_lung_cancer_poa = TRUE)

# ---------- Lung resection procedures ----------
lung_resection <- patient_proc %>%
  filter(!is.na(proc_time)) %>%
  filter(procedure_code %in% lung_proc_codes) %>%
  transmute(
    hospitalization_id,
    proc_time,
    procedure_code,
    procedure_code_format
  )

# ---------- ICU admissions (all segments) ----------
icu_adm <- icu_segments %>%
  transmute(
    hospitalization_id,
    icu_in_time  = in_dttm,
    icu_out_time = out_dttm
  ) %>%
  filter(!is.na(icu_in_time))

# ---------- Link procedure to ICU (procedure before ICU in) ----------
# Create all proc x ICU pairs per hospitalization, then keep the closest ICU AFTER proc.
proc_icu_pairs <- lung_resection %>%
  inner_join(icu_adm, by = "hospitalization_id") %>%
  mutate(
    proc_to_icu_h = as.numeric(difftime(icu_in_time, proc_time, units = "hours"))
  ) %>%
  filter(proc_to_icu_h >= 0) %>%
  filter(is.infinite(PROC_TO_ICU_MAX_H) | proc_to_icu_h <= PROC_TO_ICU_MAX_H) %>%
  group_by(hospitalization_id) %>%
  slice_min(proc_to_icu_h, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---------- Cohort builder (single cohort) ----------
build_lung_resection_icu <- function(base_df) {
  
  df <- base_df %>%
    left_join(lung_poa, by = "hospitalization_id") %>%
    left_join(proc_icu_pairs, by = "hospitalization_id") %>%
    mutate(
      has_lung_cancer_poa = coalesce(has_lung_cancer_poa, FALSE),
      has_resection_icu_link = !is.na(proc_time) & !is.na(icu_in_time),
      include =
        coalesce(adult, FALSE) &
        coalesce(has_demo, FALSE) &
        coalesce(has_geo, FALSE) &
        coalesce(has_lung_cancer_poa, FALSE) &
        coalesce(has_resection_icu_link, FALSE)
    )
  
  # Set t0
  df <- df %>%
    mutate(
      t0 = case_when(
        T0_ANCHOR == "icu_in" ~ icu_in_time,
        T0_ANCHOR == "proc"   ~ proc_time,
        TRUE ~ icu_in_time
      )
    )
  
  cohort <- df %>%
    filter(include) %>%
    transmute(
      cohort = "lung_resection_icu",
      patient_id, hospitalization_id,
      admission_dttm, discharge_dttm,
      first_icu_in, last_icu_out, icu_los_hours,
      age_years, sex_category, race_category, ethnicity_category,
      census_tract, county_code, zipcode_five_digit, zipcode_nine_digit,
      proc_time, icu_in_time, icu_out_time,
      proc_to_icu_h,
      procedure_code, procedure_code_format,
      t0
    )
  
  exclusions <- df %>%
    filter(!include) %>%
    mutate(reason = case_when(
      !coalesce(adult, FALSE) ~ "Under 18 or missing age",
      !coalesce(has_demo, FALSE) ~ "Missing demographics",
      !coalesce(has_geo, FALSE)  ~ "Missing geo code",
      !coalesce(has_lung_cancer_poa, FALSE) ~ "No lung cancer dx present-on-admission (poa_present != 1 or no matching dx code)",
      !coalesce(has_resection_icu_link, FALSE) ~ glue("No qualifying resection before ICU (proc<icu) within {PROC_TO_ICU_MAX_H}h window"),
      TRUE ~ "Other"
    )) %>%
    select(patient_id, hospitalization_id, reason)
  
  # Flow table (reviewer-friendly and debugging-friendly)
  step0 <- df
  step1 <- step0 %>% filter(coalesce(adult, FALSE))
  step2 <- step1 %>% filter(coalesce(has_demo, FALSE))
  step3 <- step2 %>% filter(coalesce(has_geo, FALSE))
  step4 <- step3 %>% filter(coalesce(has_lung_cancer_poa, FALSE))
  step5 <- step4 %>% filter(coalesce(has_resection_icu_link, FALSE))
  
  flow <- tibble(
    step = c(
      "ICU candidates (date range)",
      glue(">= {ADULT_AGE_YEARS} years"),
      "Demographics present",
      "Geography present",
      "Lung cancer POA present",
      glue("Resection before ICU (<= {PROC_TO_ICU_MAX_H}h after procedure)")
    ),
    remaining = c(nrow(step0), nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5))
  ) %>%
    mutate(excluded_at_step = lag(remaining, default = remaining[1]) - remaining)
  
  list(cohort = cohort, exclusions = exclusions, flow = flow)
}

res <- build_lung_resection_icu(base)

cohort_lung <- res$cohort
excluded_lung <- res$exclusions
flow_lung <- res$flow

cat("\nCohort selection summary:\n")
cat("  ICU candidates:                ", nrow(base), "\n", sep = "")
cat("  Lung resection ICU included:   ", nrow(cohort_lung), "\n", sep = "")
cat("  Excluded:                      ", nrow(excluded_lung), "\n", sep = "")

rs_raw <- clif_tables[["clif_respiratory_support"]] %>%
  rename_with(tolower) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    device_category = as.character(device_category)
  )

# Join to your lung resection ICU cohort
rs_at_t0 <- rs_raw %>%
  inner_join(
    cohort_lung %>%
      transmute(hospitalization_id,
                t0 = icu_in_time),
    by = "hospitalization_id"
  ) %>%
  mutate(
    dt_hours = as.numeric(difftime(recorded_dttm, t0, units = "hours"))
  ) %>%
  filter(abs(dt_hours) <= 2) %>%  # within ±2h of ICU admission
  group_by(hospitalization_id) %>%
  slice_min(abs(dt_hours), n = 1, with_ties = FALSE) %>%
  ungroup()

rs_at_t0 <- rs_at_t0 %>%
  mutate(
    rs_group = case_when(
      device_category == "IMV" ~ "Invasive Mechanical Ventilation",
      device_category %in% c("NIPPV","CPAP") ~ "Noninvasive Ventilation",
      device_category == "High Flow NC" ~ "High Flow Nasal Cannula",
      is.na(device_category) ~ "No recorded support",
      TRUE ~ "Other"
    )
  )

rs_summary <- rs_at_t0 %>%
  count(rs_group) %>%
  mutate(percent = n / sum(n) * 100)

rs_pre_icu <- rs_raw %>%
  inner_join(
    cohort_lung %>%
      select(hospitalization_id, proc_time, icu_in_time),
    by = "hospitalization_id"
  ) %>%
  filter(device_category == "IMV",
         recorded_dttm < icu_in_time,
         recorded_dttm > proc_time - dhours(12)) %>%
  distinct(hospitalization_id) %>%
  mutate(imv_pre_icu = TRUE)

# ==========================
# Respiratory support visuals
# ==========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(readr)
})

# ---------- Helper: consistent theme ----------
theme_pub <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(p, filename, width = 9, height = 5.5, dpi = 320) {
  path <- file.path(out_dir, filename)
  ggsave(path, plot = p, width = width, height = height, dpi = dpi, units = "in")
  message("Saved: ", path)
}

# ---------- Define respiratory support group mapping ----------
map_rs_group <- function(device_category) {
  d <- as.character(device_category)
  case_when(
    d == "IMV" ~ "IMV",
    d %in% c("NIPPV","CPAP") ~ "NIV/CPAP",
    d == "High Flow NC" ~ "HFNC",
    d %in% c("NC","Nasal Cannula","Simple Mask","Non-Rebreather","Face Mask","Venturi") ~ "Low-flow O2",
    d %in% c("Room Air","RA") ~ "Room Air",
    is.na(d) | d == "" ~ "No RS recorded",
    TRUE ~ "Other"
  )
}

# ---------- 0) Build rs_at_t0 if not already present ----------
# Requires: rs_raw (respiratory_support) and cohort_lung (with icu_in_time)
# This is safe to rerun; it will overwrite rs_at_t0.
rs_at_t0 <- rs_raw %>%
  inner_join(
    cohort_lung %>% transmute(hospitalization_id, t0 = icu_in_time),
    by = "hospitalization_id"
  ) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    dt_hours = as.numeric(difftime(recorded_dttm, t0, units = "hours"))
  ) %>%
  filter(abs(dt_hours) <= 2) %>%                             # ±2h window around ICU entry
  group_by(hospitalization_id) %>%
  slice_min(abs(dt_hours), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    rs_group = map_rs_group(device_category),
    fio2 = suppressWarnings(as.numeric(fio2_set))
  )

# Ensure a stable ordering for plots
rs_levels <- c("Room Air","Low-flow O2","HFNC","NIV/CPAP","IMV","Other","No RS recorded")
rs_at_t0 <- rs_at_t0 %>%
  mutate(rs_group = factor(rs_group, levels = rs_levels))

# ---------- 1) Bar chart: device category at ICU entry ----------
rs_summary <- rs_at_t0 %>%
  count(rs_group, name = "n") %>%
  mutate(pct = n / sum(n))

p1 <- ggplot(rs_summary, aes(x = rs_group, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(n, "\n", percent(pct, accuracy = 0.1))),
            vjust = -0.2, size = 3.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Respiratory support at ICU admission (t0)",
    x = NULL, y = "Hospitalizations"
  ) +
  theme_pub(14) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p1, "fig_rs_at_icu_entry_bar.png", width = 9.5, height = 5.5)

# Also export the table used for the plot
readr::write_csv(rs_summary, file.path(out_dir, "table_rs_at_icu_entry.csv"))

# ---------- 2) By exposure quantile (if exposure exists in cohort_lung) ----------
# Expected column names (edit to your actual exposure columns):
#   pm25_1y, no2_1y, pm25, no2, etc.
# exposure_candidates <- c("pm25_1y","no2_1y","pm25","no2","pm2_5","no2_ppb")
# expo_name <- exposure_candidates[exposure_candidates %in% names(cohort_lung)][1]
# 
# if (!is.na(expo_name)) {
#   
#   dat_expo <- rs_at_t0 %>%
#     inner_join(cohort_lung %>% select(hospitalization_id, all_of(expo_name)),
#                by = "hospitalization_id") %>%
#     mutate(
#       exposure = suppressWarnings(as.numeric(.data[[expo_name]])),
#       expo_q = ntile(exposure, 4),
#       expo_q = factor(expo_q, levels = 1:4,
#                       labels = c("Q1 (lowest)","Q2","Q3","Q4 (highest)"))
#     ) %>%
#     filter(!is.na(exposure))
#   
#   rs_by_q <- dat_expo %>%
#     count(expo_q, rs_group) %>%
#     group_by(expo_q) %>%
#     mutate(prop = n / sum(n)) %>%
#     ungroup()
#   
#   p2 <- ggplot(rs_by_q, aes(x = expo_q, y = prop, fill = rs_group)) +
#     geom_col(position = "fill") +
#     scale_y_continuous(labels = percent_format()) +
#     labs(
#       title = glue("Respiratory support at ICU admission by {expo_name} quartile"),
#       x = NULL, y = "Proportion", fill = "Resp support"
#     ) +
#     theme_pub(14)
#   
#   save_plot(p2, glue("fig_rs_by_{expo_name}_quartile_stacked.png"), width = 9.5, height = 5.8)
#   readr::write_csv(rs_by_q, file.path(out_dir, glue("table_rs_by_{expo_name}_quartile.csv")))
#   
# } else {
#   message("No exposure column found in cohort_lung for stratified plot. Skipping p2.")
# }

# ---------- 3) Heatmap: rs_group x FiO2 at t0 (if FiO2 available) ----------
# FiO2 often stored as 0-1 or 21-100; we’ll standardize to percent.
rs_fio2 <- rs_at_t0 %>%
  mutate(
    fio2_pct = case_when(
      is.na(fio2) ~ NA_real_,
      fio2 <= 1 ~ fio2 * 100,
      fio2 > 1 ~ fio2
    ),
    fio2_bin = cut(fio2_pct,
                   breaks = c(-Inf, 30, 40, 50, 60, 70, 80, 90, Inf),
                   right = TRUE,
                   labels = c("≤30","31–40","41–50","51–60","61–70","71–80","81–90",">90"))
  ) %>%
  filter(!is.na(fio2_bin))

if (nrow(rs_fio2) > 0) {
  hm <- rs_fio2 %>%
    count(rs_group, fio2_bin) %>%
    group_by(rs_group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  p3 <- ggplot(hm, aes(x = fio2_bin, y = rs_group, fill = prop)) +
    geom_tile() +
    scale_fill_continuous(labels = percent_format()) +
    labs(
      title = "Initial FiO₂ distribution at ICU admission by respiratory support",
      x = "FiO₂ at/near t0 (binned)", y = NULL, fill = "Row %"
    ) +
    theme_pub(14)
  
  save_plot(p3, "fig_rs_x_fio2_heatmap.png", width = 9.5, height = 5.8)
  readr::write_csv(hm, file.path(out_dir, "table_rs_x_fio2_heatmap.csv"))
} else {
  message("FiO2 not available in rs_at_t0 (fio2_set missing/empty). Skipping p3.")
}

# ---------- 4) Escalation within 24h: device at t0 -> max device within 24h ----------
# Define an ordinal severity for max device determination
device_severity <- function(device_category) {
  d <- as.character(device_category)
  case_when(
    d == "IMV" ~ 5L,
    d %in% c("NIPPV","CPAP") ~ 4L,
    d == "High Flow NC" ~ 3L,
    d %in% c("NC","Nasal Cannula","Simple Mask","Non-Rebreather","Face Mask","Venturi") ~ 2L,
    d %in% c("Room Air","RA") ~ 1L,
    TRUE ~ 0L
  )
}

# Build t0 per hospitalization (ICU in time)
t0_tbl <- cohort_lung %>% transmute(hospitalization_id, t0 = icu_in_time)

rs_24h <- rs_raw %>%
  inner_join(t0_tbl, by = "hospitalization_id") %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    dt_h = as.numeric(difftime(recorded_dttm, t0, units = "hours"))
  ) %>%
  filter(dt_h >= 0, dt_h <= 24) %>%
  mutate(
    rs_group = map_rs_group(device_category),
    sev = device_severity(device_category)
  )

# Determine max support within 24h by severity (ties resolved by earliest time at that severity)
rs_max24 <- rs_24h %>%
  group_by(hospitalization_id) %>%
  filter(sev == max(sev, na.rm = TRUE)) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(hospitalization_id, rs_max24 = factor(rs_group, levels = rs_levels))

# Join baseline (t0) group
transitions <- rs_at_t0 %>%
  select(hospitalization_id, rs_t0 = rs_group) %>%
  left_join(rs_max24, by = "hospitalization_id") %>%
  mutate(
    rs_max24 = fct_na_value_to_level(rs_max24, level = "No RS recorded"),
    rs_t0    = fct_na_value_to_level(rs_t0, level = "No RS recorded")
  )

# Transition matrix heatmap (nice and grant-friendly)
mat <- transitions %>%
  count(rs_t0, rs_max24) %>%
  group_by(rs_t0) %>%
  mutate(row_prop = n / sum(n)) %>%
  ungroup()

p4 <- ggplot(mat, aes(x = rs_max24, y = rs_t0, fill = row_prop)) +
  geom_tile() +
  geom_text(aes(label = n), size = 3.5) +
  scale_fill_continuous(labels = percent_format()) +
  labs(
    title = "Respiratory support transitions (t0 → max within 24h of ICU admission)",
    x = "Max support within 24h",
    y = "Support at ICU admission (t0)",
    fill = "Row %"
  ) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p4, "fig_rs_transition_t0_to_max24_heatmap.png", width = 10.5, height = 6.2)
readr::write_csv(mat, file.path(out_dir, "table_rs_transition_t0_to_max24.csv"))

message("Respiratory support visualizations complete.")

      



















# # ---------- Save outputs ----------
# SITE_NAME   <- sanitize_tag(site_name)
# SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")
# 
# make_name <- function(result_name, ext = "csv") {
#   paste0(result_name, "_", SITE_NAME, "_", SYSTEM_DATE, ".", ext)
# }
# 
# out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE))
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# write_csv(cohort_lung,   file.path(out_dir, make_name("cohort_lung_resection_icu")))
# write_csv(excluded_lung, file.path(out_dir, make_name("exclusion_lung_resection_icu")))
# write_csv(flow_lung,     file.path(out_dir, make_name("flow_lung_resection_icu")))
# 
# cat("\nSaved outputs to: ", out_dir, "\n", sep = "")
# cat("  - ", make_name("cohort_lung_resection_icu"), "\n", sep = "")
# cat("  - ", make_name("exclusion_lung_resection_icu"), "\n", sep = "")
# cat("  - ", make_name("flow_lung_resection_icu"), "\n", sep = "")
# 
# # ---------- Optional: keep minimal objects ----------
# keep_vars <- c("clif_tables", "cohort_lung", "repo", "out_dir", "flow_lung", "excluded_lung")
# rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)
# 
# message("\nCohort identification complete.")
