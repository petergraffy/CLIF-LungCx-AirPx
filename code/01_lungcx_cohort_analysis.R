
# ================================================================================================
# CLIF | ICU Post–Lung Cancer Resection Cohort + Exposome + Respiratory Support Summary Analyses
# PI: Peter Graffy
#
# Goal (federated-friendly):
#   1) Build cohort (dx POA + curated resection codes + ICU after procedure)
#   2) Link county-year PM2.5 / NO2 exposures (1y + 3y pre-op)
#   3) Derive respiratory support course features and IMV hours (ICU + first 72h)
#   4) Produce small-sample-appropriate summaries + exact tests + effect size tables
#   5) Save ALL outputs into repo/output/run_[SITE]_[DATE]/ for cross-site pooling
#
# Notes:
#   - Avoids relying on procedure_code_format for code system inference.
#   - Uses nonparametric tests and exact tests (n~56).
#   - Produces per-site tables + metadata for federated meta-analysis.
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(lubridate)
  library(readr)
  library(glue)
  library(fst)
  library(arrow)
  library(scales)
})

# --------------------------- Config ---------------------------
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

# --------------------------- Parameters ---------------------------
START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")

ADULT_AGE_YEARS   <- 18
PROC_TO_ICU_MAX_H <- 48
T0_ANCHOR         <- "icu_in"  # "icu_in" or "proc"

# dx prefixes live here (ICD9/10-CM lung cancer prefixes)
CODES_PATH <- file.path(repo, "outlier-thresholds", "table.tsv")
stopifnot(file.exists(CODES_PATH))

# exposome input files (county-year)
EXPO_DIR <- file.path(repo, "exposome")
PM25_PATH <- file.path(EXPO_DIR, "pm25_county_year.csv")
NO2_PATH  <- file.path(EXPO_DIR, "no2_county_year.csv")
stopifnot(file.exists(PM25_PATH), file.exists(NO2_PATH))

# --------------------------- Helpers ---------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

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
norm_code <- function(x) {
  x <- toupper(as.character(x))
  str_replace_all(x, "[^A-Z0-9]", "")
}
code_matches_any_prefix <- function(code_vec, prefixes) {
  code_vec <- norm_code(code_vec)
  prefixes <- norm_code(prefixes)
  vapply(code_vec, function(cd) any(str_starts(cd, prefixes)), logical(1))
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
norm_county <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^0-9]", "")
  x <- ifelse(nchar(x) == 4, paste0("0", x), x) # pad to 5 if needed
  x
}

theme_pub <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# --------------------------- Output folder + I/O ---------------------------
SITE_NAME   <- sanitize_tag(site_name)
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")

out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

make_name <- function(stub, ext="csv") paste0(stub, "_", SITE_NAME, "_", SYSTEM_DATE, ".", ext)

save_csv <- function(df, stub) {
  p <- file.path(out_dir, make_name(stub, "csv"))
  readr::write_csv(df, p)
  message("Saved: ", p)
}
save_rds <- function(obj, stub) {
  p <- file.path(out_dir, make_name(stub, "rds"))
  saveRDS(obj, p)
  message("Saved: ", p)
}
save_plot <- function(p, stub, w=9, h=5.5, dpi=320) {
  fp <- file.path(out_dir, make_name(stub, "png"))
  ggsave(fp, plot = p, width = w, height = h, dpi = dpi, units = "in")
  message("Saved: ", fp)
}

# Federated analysis metadata (site-level)
analysis_meta <- tibble(
  site_name = SITE_NAME,
  run_date  = SYSTEM_DATE,
  start_date = as.character(START_DATE),
  end_date   = as.character(END_DATE),
  adult_age_years = ADULT_AGE_YEARS,
  proc_to_icu_max_h = PROC_TO_ICU_MAX_H,
  t0_anchor = T0_ANCHOR
)
save_csv(analysis_meta, "meta_run_parameters")

# --------------------------- Load CLIF tables ---------------------------
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
if (length(all_files) == 0) stop("No CLIF files found under: ", tables_path)

bn <- basename(all_files)
looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)
base_no_ext <- tools::file_path_sans_ext(tolower(bn))
base_norm <- ifelse(looks_clif, base_no_ext, paste0("clif_", base_no_ext))
found_map <- stats::setNames(all_files, base_norm)

required_raw <- c("patient","hospitalization","adt","hospital_diagnosis","patient_procedures","respiratory_support","labs","vitals")
required_files <- paste0("clif_", required_raw)
missing <- setdiff(required_files, names(found_map))
if (length(missing) > 0) {
  cat("Detected CLIF-like files:\n"); print(sort(unique(names(found_map))))
  stop("Missing required tables: ", paste(missing, collapse = ", "))
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext, " for ", path))
}

clif_paths  <- found_map[required_files]
clif_tables <- lapply(clif_paths, read_any)
names(clif_tables) <- required_files
cat("Loaded tables: ", paste(names(clif_tables), collapse = ", "), "\n")

get_min <- function(tbl_name, cols) {
  nm <- paste0("clif_", tbl_name)
  out <- clif_tables[[nm]] %>% rename_with(tolower)
  cols_keep <- intersect(tolower(cols), names(out))
  out %>% select(any_of(cols_keep))
}

# --------------------------- Load dx prefixes ---------------------------
codes <- readr::read_tsv(CODES_PATH, show_col_types = FALSE) %>%
  rename_with(tolower)
stopifnot(all(c("code type","code") %in% names(codes)))

lung_dx_prefixes <- codes %>%
  mutate(
    code_system = toupper(str_replace_all(`code type`, "[^A-Z0-9]", "")),
    code = norm_code(code)
  ) %>%
  filter(code_system %in% c("ICD10CM","ICD9CM")) %>%
  pull(code) %>%
  unique()

save_csv(tibble(lung_dx_prefix=lung_dx_prefixes), "qc_lung_dx_prefixes")

# --------------------------- Curated resection procedure codes ---------------------------
resection_cpt <- norm_code(c("32480","32482","32663","32484","32505","32666","32669","32440","32442"))
resection_icd10pcs <- norm_code(c(
  "0BTC0ZZ","0BTC4ZZ","0BTD0ZZ","0BTD4ZZ","0BTG0ZZ","0BTG4ZZ",
  "0BTH0ZZ","0BTH4ZZ","0BTJ0ZZ","0BTJ4ZZ","0BTK0ZZ","0BTK4ZZ"
))

surg_from_code <- function(proc_code) {
  pc <- norm_code(proc_code)
  case_when(
    pc %in% c("32440","32442") ~ "Pneumonectomy",
    pc %in% c("32480","32482","32663") ~ "Lobectomy",
    pc %in% c("32484","32669") ~ "Segmentectomy",
    pc %in% c("32505","32666") ~ "Wedge",
    str_starts(pc, "0BTJ") ~ "Pneumonectomy",
    str_starts(pc, "0BTK") ~ "Pneumonectomy",
    str_starts(pc, "0BTC") ~ "Lobectomy",
    str_starts(pc, "0BTD") ~ "Lobectomy",
    str_starts(pc, "0BTG") ~ "Lobectomy",
    str_starts(pc, "0BTH") ~ "Lobectomy",
    TRUE ~ "Other"
  )
}

save_csv(
  tibble(code=c(resection_cpt, resection_icd10pcs),
         code_type=c(rep("CPT", length(resection_cpt)), rep("ICD10PCS", length(resection_icd10pcs))),
         surgery_category=surg_from_code(code)),
  "qc_resection_code_list"
)

# --------------------------- Minimal tables ---------------------------
patient <- get_min("patient", c("patient_id","birth_date","sex_category","race_category","ethnicity_category")) %>%
  mutate(birth_date = safe_date(birth_date))

hospitalization <- get_min("hospitalization", c(
  "patient_id","hospitalization_id","admission_dttm","discharge_dttm","discharge_category",
  "age_at_admission","zipcode_nine_digit","zipcode_five_digit","census_tract","county_code"
)) %>%
  mutate(
    admission_dttm = safe_posix(admission_dttm),
    discharge_dttm = safe_posix(discharge_dttm),
    discharge_category = as.character(discharge_category)
  ) %>%
  mutate(
    discharge_l = tolower(coalesce(discharge_category, "")),
    death_in_hosp = str_detect(discharge_l, "expired|death"),
    hospice_discharge = str_detect(discharge_l, "hospice"),
    death_or_hospice = death_in_hosp | hospice_discharge
  ) %>%
  select(-discharge_l)

adt <- get_min("adt", c("hospitalization_id","in_dttm","out_dttm","location_category","location_type")) %>%
  mutate(
    in_dttm  = safe_posix(in_dttm),
    out_dttm = safe_posix(out_dttm),
    location_category = as.character(location_category)
  )

hospital_dx <- get_min("hospital_diagnosis", c(
  "hospitalization_id","diagnosis_code","diagnosis_code_format","poa_present","diagnosis_primary"
)) %>%
  mutate(
    icd_code = norm_code(diagnosis_code),
    poa_present = suppressWarnings(as.integer(poa_present))
  )

patient_proc <- get_min("patient_procedures", c(
  "hospitalization_id","procedure_code","procedure_code_format","procedure_billed_dttm"
)) %>%
  mutate(
    procedure_code = norm_code(procedure_code),
    procedure_billed_dttm = safe_posix(procedure_billed_dttm),
    proc_time = procedure_billed_dttm
  )

rs_raw <- clif_tables[["clif_respiratory_support"]] %>%
  rename_with(tolower) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    device_category = as.character(device_category)
  )

# --------------------------- Cohort Build ---------------------------
# A) Base hospitalizations in window + adults + geo
base <- hospitalization %>%
  filter(!is.na(admission_dttm),
         admission_dttm >= START_DATE,
         admission_dttm <= END_DATE) %>%
  left_join(patient %>% select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
            by = "patient_id") %>%
  mutate(
    age_years = coalesce(
      suppressWarnings(as.numeric(age_at_admission)),
      ifelse(!is.na(birth_date),
             as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
    ),
    has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
    adult    = !is.na(age_years) & age_years >= ADULT_AGE_YEARS,
    has_geo  = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code)
  )

# B) Lung cancer dx POA=1
lung_poa <- hospital_dx %>%
  filter(poa_present == 1) %>%
  filter(code_matches_any_prefix(icd_code, lung_dx_prefixes)) %>%
  distinct(hospitalization_id) %>%
  mutate(has_lung_cancer_poa = TRUE)

# C) Resection procedure match by code shape (not procedure_code_format)
proc_resection <- patient_proc %>%
  filter(!is.na(proc_time), !is.na(procedure_code), procedure_code != "") %>%
  mutate(
    is_cpt = str_detect(procedure_code, "^[0-9]+$"),
    is_pcs = str_detect(procedure_code, "^[A-Z0-9]+$") & str_starts(procedure_code, "0B")
  ) %>%
  filter(
    (is_cpt & procedure_code %in% resection_cpt) |
      (is_pcs & procedure_code %in% resection_icd10pcs)
  ) %>%
  mutate(surgery_category = surg_from_code(procedure_code)) %>%
  transmute(hospitalization_id, proc_time, procedure_code, surgery_category)

# D) ICU segments
icu_segments <- adt %>%
  mutate(is_icu = str_detect(tolower(coalesce(location_category, "")), "icu")) %>%
  filter(is_icu, !is.na(in_dttm)) %>%
  transmute(hospitalization_id,
            icu_in_time = in_dttm,
            icu_out_time = out_dttm)

icu_bounds <- icu_segments %>%
  group_by(hospitalization_id) %>%
  summarize(
    first_icu_in = suppressWarnings(min(icu_in_time, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(icu_out_time, na.rm = TRUE)),
    .groups="drop"
  ) %>%
  mutate(
    first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
    last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units="hours"))
  )

# E) eligible (dx + resection)
eligible <- base %>%
  inner_join(lung_poa, by="hospitalization_id") %>%
  inner_join(proc_resection, by="hospitalization_id")

# F) link proc -> closest ICU admission after proc within window
proc_icu_pairs <- eligible %>%
  inner_join(icu_segments, by="hospitalization_id") %>%
  mutate(proc_to_icu_h = as.numeric(difftime(icu_in_time, proc_time, units="hours"))) %>%
  filter(proc_to_icu_h >= 0) %>%
  filter(is.infinite(PROC_TO_ICU_MAX_H) | proc_to_icu_h <= PROC_TO_ICU_MAX_H) %>%
  group_by(hospitalization_id) %>%
  slice_min(proc_to_icu_h, n=1, with_ties = FALSE) %>%
  ungroup()

cohort_lung <- proc_icu_pairs %>%
  left_join(icu_bounds, by="hospitalization_id") %>%
  mutate(
    t0 = case_when(
      T0_ANCHOR == "icu_in" ~ icu_in_time,
      T0_ANCHOR == "proc"   ~ proc_time,
      TRUE ~ icu_in_time
    )
  ) %>%
  transmute(
    cohort = "lung_resection_icu",
    patient_id, hospitalization_id,
    admission_dttm, discharge_dttm,
    discharge_category, death_in_hosp, hospice_discharge, death_or_hospice,
    first_icu_in, last_icu_out, icu_los_hours,
    age_years, sex_category, race_category, ethnicity_category,
    census_tract, county_code, zipcode_five_digit, zipcode_nine_digit,
    proc_time, procedure_code, surgery_category,
    icu_in_time, icu_out_time, proc_to_icu_h,
    t0,
    adult, has_demo, has_geo
  ) %>%
  mutate(include = coalesce(adult, FALSE) & coalesce(has_demo, FALSE) & coalesce(has_geo, FALSE)) %>%
  filter(include) %>%
  select(-include)

# Flow + exclusions (federated-friendly)
step0 <- base
step1 <- step0 %>% filter(adult)
step2 <- step1 %>% filter(has_demo)
step3 <- step2 %>% filter(has_geo)
step4 <- step3 %>% inner_join(lung_poa, by="hospitalization_id")
step5 <- step4 %>% inner_join(proc_resection %>% distinct(hospitalization_id), by="hospitalization_id")
step6 <- step5 %>% inner_join(icu_segments %>% distinct(hospitalization_id), by="hospitalization_id")

flow_lung <- tibble(
  step = c(
    "Hospitalizations in date window",
    glue(">= {ADULT_AGE_YEARS} years"),
    "Demographics present",
    "Geography present",
    "Lung cancer POA present",
    "Resection procedure present",
    "Any ICU segment present"
  ),
  remaining = c(nrow(step0), nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5), nrow(step6))
) %>%
  mutate(excluded_at_step = lag(remaining, default = remaining[1]) - remaining)

exclusion_lung <- base %>%
  mutate(
    reason = case_when(
      !adult ~ "Under 18 or missing age",
      !has_demo ~ "Missing demographics",
      !has_geo ~ "Missing geo code",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(lung_poa %>% transmute(hospitalization_id, dx_ok=TRUE), by="hospitalization_id") %>%
  left_join(proc_resection %>% distinct(hospitalization_id) %>% mutate(proc_ok=TRUE), by="hospitalization_id") %>%
  left_join(icu_segments %>% distinct(hospitalization_id) %>% mutate(icu_ok=TRUE), by="hospitalization_id") %>%
  mutate(
    reason = coalesce(
      reason,
      ifelse(is.na(dx_ok), "No lung cancer dx POA=1 matching prefixes", NA_character_),
      ifelse(is.na(proc_ok), "No qualifying lung resection procedure code", NA_character_),
      ifelse(is.na(icu_ok), "No ICU stay (no ADT segment with ICU location_category)", NA_character_),
      "Other"
    )
  ) %>%
  filter(!(hospitalization_id %in% cohort_lung$hospitalization_id)) %>%
  select(patient_id, hospitalization_id, reason)

cat("\nCohort selection summary:\n")
cat("  Base hospitalizations:         ", nrow(step0), "\n", sep="")
cat("  Eligible (dx+proc):            ", nrow(step5), "\n", sep="")
cat("  Eligible + ICU:                ", nrow(step6), "\n", sep="")
cat("  Final cohort included:         ", nrow(cohort_lung), "\n", sep="")

# --------------------------- Exposome linkage (county-year) ---------------------------
pm25 <- readr::read_csv(PM25_PATH, show_col_types = FALSE)
no2  <- readr::read_csv(NO2_PATH,  show_col_types = FALSE)

# IMPORTANT: adjust these renames to your actual columns ONCE, then keep stable across sites.
# Below assumes (GEOID, year, *_mean). If not, change here and nowhere else.
pm25 <- pm25 %>%
  rename(county = GEOID, year = year, pm25 = pm25_mean) %>%
  mutate(county = norm_county(county), year = as.integer(year), pm25 = as_num(pm25))

no2 <- no2 %>%
  rename(county = GEOID, year = year, no2 = no2_mean) %>%
  mutate(county = norm_county(county), year = as.integer(year), no2 = as_num(no2))

cohort_expo <- cohort_lung %>%
  mutate(
    county = norm_county(county_code),
    proc_date = as.Date(proc_time),
    proc_year = as.integer(format(proc_date, "%Y")),
    y1 = proc_year - 1,
    y3_start = proc_year - 3,
    y3_end   = proc_year - 1
  )

pm25_1y <- cohort_expo %>%
  select(hospitalization_id, county, y1) %>%
  left_join(pm25 %>% select(county, year, pm25), by = c("county"="county", "y1"="year")) %>%
  rename(pm25_1y = pm25) %>%
  select(hospitalization_id, pm25_1y)

no2_1y <- cohort_expo %>%
  select(hospitalization_id, county, y1) %>%
  left_join(no2 %>% select(county, year, no2), by = c("county"="county", "y1"="year")) %>%
  rename(no2_1y = no2) %>%
  select(hospitalization_id, no2_1y)

pm25_3y <- cohort_expo %>%
  select(hospitalization_id, county, y3_start, y3_end) %>%
  left_join(pm25 %>% select(county, year, pm25), by="county") %>%
  filter(year >= y3_start, year <= y3_end) %>%
  group_by(hospitalization_id) %>%
  summarise(pm25_3y = mean(pm25, na.rm=TRUE), .groups="drop")

no2_3y <- cohort_expo %>%
  select(hospitalization_id, county, y3_start, y3_end) %>%
  left_join(no2 %>% select(county, year, no2), by="county") %>%
  filter(year >= y3_start, year <= y3_end) %>%
  group_by(hospitalization_id) %>%
  summarise(no2_3y = mean(no2, na.rm=TRUE), .groups="drop")

dat <- cohort_lung %>%
  left_join(pm25_1y, by="hospitalization_id") %>%
  left_join(no2_1y,  by="hospitalization_id") %>%
  left_join(pm25_3y, by="hospitalization_id") %>%
  left_join(no2_3y,  by="hospitalization_id") %>%
  mutate(
    pm25_tertile = ifelse(!is.na(pm25_1y), ntile(pm25_1y, 3), NA_integer_),
    no2_tertile  = ifelse(!is.na(no2_1y),  ntile(no2_1y,  3), NA_integer_)
  )

qc_expo <- dat %>%
  summarise(
    n = n(),
    pm25_1y_missing = mean(is.na(pm25_1y)),
    no2_1y_missing  = mean(is.na(no2_1y)),
    pm25_3y_missing = mean(is.na(pm25_3y)),
    no2_3y_missing  = mean(is.na(no2_3y))
  )

save_csv(qc_expo, "qc_exposure_linkage")

p_expo <- dat %>%
  pivot_longer(c(pm25_1y, no2_1y), names_to="exposure", values_to="value") %>%
  ggplot(aes(x=exposure, y=value)) +
  geom_boxplot() +
  labs(title="Pre-op county-level exposures (1y)", x=NULL, y=NULL) +
  theme_pub()
save_plot(p_expo, "fig_exposure_boxplot", w=7, h=5)

# --------------------------- Respiratory support panel (federated-ready) ---------------------------
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
severity_state <- function(state) {
  case_when(
    state == "IMV" ~ 5L,
    state == "NIV/CPAP" ~ 4L,
    state == "HFNC" ~ 3L,
    state == "Low-flow O2" ~ 2L,
    state == "Room Air" ~ 1L,
    TRUE ~ 0L
  )
}

t0_tbl <- dat %>% transmute(hospitalization_id, t0 = icu_in_time)

# Panel within -2..72h around ICU in (t0)
rs_panel_72h <- rs_raw %>%
  inner_join(t0_tbl, by="hospitalization_id") %>%
  filter(!is.na(recorded_dttm), !is.na(t0)) %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    state = map_rs_group(device_category)
  ) %>%
  filter(dt_h >= -2, dt_h <= 72) %>%
  arrange(hospitalization_id, recorded_dttm)

rs_panel_export <- rs_panel_72h %>%
  arrange(hospitalization_id, dt_h) %>%
  group_by(hospitalization_id) %>%
  mutate(
    gen_patient_id = cur_group_id()
  ) %>%
  ungroup()

rs_panel_export <- rs_panel_export %>%
  select(
    gen_patient_id,
    dt_h,
    state,
    device_category
  )

save_csv(
  rs_panel_export,
  "panel_respiratory_support_-2to72h_federated"
)

# --------------------------- Course features (small n, high yield) ---------------------------
course_features <- rs_panel_72h %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_t0 = any(state == "IMV" & abs(dt_h) <= 2, na.rm=TRUE),
    imv_any_72 = any(state == "IMV" & dt_h >= 0 & dt_h <= 72, na.rm=TRUE),
    
    # time to first IMV (if any)
    imv_start_h = if (any(state=="IMV" & dt_h>=0 & dt_h<=72, na.rm=TRUE)) {
      min(dt_h[state=="IMV" & dt_h>=0 & dt_h<=72], na.rm=TRUE)
    } else NA_real_,
    
    # extubation proxy: first non-IMV after first IMV (post t0)
    extub_h = {
      idx_imv <- which(state=="IMV" & dt_h>=0)
      if (length(idx_imv)==0) NA_real_ else {
        first_imv <- idx_imv[1]
        idx_post <- which(seq_along(state) > first_imv & dt_h>=0 & state != "IMV")
        if (length(idx_post)==0) NA_real_ else dt_h[idx_post[1]]
      }
    },
    extub_24h = !is.na(extub_h) & extub_h <= 24,
    extub_48h = !is.na(extub_h) & extub_h <= 48,
    
    reintub_48h = {
      if (is.na(extub_h)) FALSE else {
        any(state=="IMV" & dt_h > extub_h & dt_h <= extub_h + 48, na.rm=TRUE)
      }
    },
    
    max_sev_24h = {
      d24 <- dt_h >= 0 & dt_h <= 24
      if (!any(d24, na.rm=TRUE)) NA_integer_ else max(severity_state(state[d24]), na.rm=TRUE)
    },
    
    .groups="drop"
  ) %>%
  mutate(
    max_support_24h = case_when(
      max_sev_24h == 5L ~ "IMV",
      max_sev_24h == 4L ~ "NIV/CPAP",
      max_sev_24h == 3L ~ "HFNC",
      max_sev_24h == 2L ~ "Low-flow O2",
      max_sev_24h == 1L ~ "Room Air",
      max_sev_24h == 0L ~ "None/Other",
      TRUE ~ NA_character_
    ),
    max_support_24h = factor(max_support_24h,
                             levels=c("Room Air","Low-flow O2","HFNC","NIV/CPAP","IMV","None/Other"))
  )

# --------------------------- IMV hours (ICU + 72h) ---------------------------
# interval approximation between successive respiratory_support records (within ICU)
rs_long <- rs_raw %>%
  inner_join(dat %>% select(hospitalization_id, icu_in_time, icu_out_time), by="hospitalization_id") %>%
  filter(!is.na(recorded_dttm), recorded_dttm >= icu_in_time, recorded_dttm <= icu_out_time) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  mutate(
    next_time = lead(recorded_dttm),
    interval_h = as.numeric(difftime(next_time, recorded_dttm, units="hours")),
    interval_h = ifelse(is.na(interval_h) | interval_h < 0, 0, interval_h)
  ) %>%
  ungroup()

imv_hours_total <- rs_long %>%
  filter(device_category == "IMV") %>%
  group_by(hospitalization_id) %>%
  summarise(imv_hours_total = sum(interval_h, na.rm=TRUE), .groups="drop")

imv_hours_72 <- rs_long %>%
  mutate(hours_from_icu = as.numeric(difftime(recorded_dttm, icu_in_time, units="hours"))) %>%
  filter(hours_from_icu >= 0, hours_from_icu <= 72) %>%
  filter(device_category == "IMV") %>%
  group_by(hospitalization_id) %>%
  summarise(imv_hours_72 = sum(interval_h, na.rm=TRUE), .groups="drop")

# Merge analytic dataset
dat2 <- dat %>%
  left_join(course_features, by="hospitalization_id") %>%
  left_join(imv_hours_total, by="hospitalization_id") %>%
  left_join(imv_hours_72, by="hospitalization_id") %>%
  mutate(
    imv_hours_total = coalesce(imv_hours_total, 0),
    imv_hours_72    = coalesce(imv_hours_72, 0),
    hosp_los_days = as.numeric(difftime(discharge_dttm, admission_dttm, units="days")),
    icu_los_days  = icu_los_hours / 24
  )

# --------------------------- Descriptive tables (federated-ready) ---------------------------
# Baseline characteristics by exposure tertile (counts + medians)
baseline_by_tertile <- dat2 %>%
  filter(!is.na(pm25_tertile)) %>%
  group_by(pm25_tertile) %>%
  summarise(
    n = n(),
    age_median = median(age_years, na.rm=TRUE),
    imv_t0_n = sum(imv_t0, na.rm=TRUE),
    imv_any72_n = sum(imv_any_72, na.rm=TRUE),
    death_n = sum(death_in_hosp, na.rm=TRUE),
    death_or_hospice_n = sum(death_or_hospice, na.rm=TRUE),
    hosp_los_median = median(hosp_los_days, na.rm=TRUE),
    icu_los_median  = median(icu_los_days, na.rm=TRUE),
    imv_hours72_median = median(imv_hours_72, na.rm=TRUE),
    imv_hours_total_median = median(imv_hours_total, na.rm=TRUE),
    .groups="drop"
  )
save_csv(baseline_by_tertile, "table_baseline_by_pm25_tertile")

# Surgery category distribution vs exposure tertile (check confounding)
surg_by_tertile <- dat2 %>%
  count(surgery_category, pm25_tertile, name="n") %>%
  group_by(surgery_category) %>%
  mutate(prop_within_surgery = n/sum(n)) %>%
  ungroup()
save_csv(surg_by_tertile, "table_surgery_by_pm25_tertile")

# --------------------------- Respiratory support differences (small n, exact/nonparametric) ---------------------------
# 1) Death vs exposure tertiles (exact)
tab_death_pm25 <- table(dat2$pm25_tertile, dat2$death_in_hosp)
tab_death_no2  <- table(dat2$no2_tertile,  dat2$death_in_hosp)

fisher_death_pm25 <- broom::tidy(fisher.test(tab_death_pm25))
fisher_death_no2  <- broom::tidy(fisher.test(tab_death_no2))

save_csv(as_tibble(tab_death_pm25, .name_repair="minimal"), "tab_death_by_pm25_tertile_raw")
save_csv(as_tibble(tab_death_no2,  .name_repair="minimal"), "tab_death_by_no2_tertile_raw")
save_csv(fisher_death_pm25, "test_fisher_death_by_pm25_tertile")
save_csv(fisher_death_no2,  "test_fisher_death_by_no2_tertile")

# 2) Trend test for death across PM2.5 tertile (interpret cautiously)
trend_death_pm25 <- prop.trend.test(
  x = tapply(dat2$death_in_hosp, dat2$pm25_tertile, sum),
  n = tapply(dat2$death_in_hosp, dat2$pm25_tertile, length)
)
trend_tbl <- tibble(
  test="prop.trend.test",
  outcome="death_in_hosp",
  exposure="pm25_tertile",
  statistic=unname(trend_death_pm25$statistic),
  df=unname(trend_death_pm25$parameter),
  p_value=trend_death_pm25$p.value
)
save_csv(trend_tbl, "test_trend_death_by_pm25_tertile")

# 3) IMV hours vs exposure tertile (Kruskal-Wallis)
kw_imv_total_pm25 <- kruskal.test(imv_hours_total ~ factor(pm25_tertile), data=dat2)
kw_imv_72_pm25    <- kruskal.test(imv_hours_72    ~ factor(pm25_tertile), data=dat2)
kw_imv_total_no2  <- kruskal.test(imv_hours_total ~ factor(no2_tertile),  data=dat2)
kw_imv_72_no2     <- kruskal.test(imv_hours_72    ~ factor(no2_tertile),  data=dat2)

kw_tbl <- tibble(
  test="kruskal.test",
  outcome=c("imv_hours_total","imv_hours_72","imv_hours_total","imv_hours_72"),
  exposure=c("pm25_tertile","pm25_tertile","no2_tertile","no2_tertile"),
  statistic=c(unname(kw_imv_total_pm25$statistic),
              unname(kw_imv_72_pm25$statistic),
              unname(kw_imv_total_no2$statistic),
              unname(kw_imv_72_no2$statistic)),
  df=c(unname(kw_imv_total_pm25$parameter),
       unname(kw_imv_72_pm25$parameter),
       unname(kw_imv_total_no2$parameter),
       unname(kw_imv_72_no2$parameter)),
  p_value=c(kw_imv_total_pm25$p.value,
            kw_imv_72_pm25$p.value,
            kw_imv_total_no2$p.value,
            kw_imv_72_no2$p.value)
)
save_csv(kw_tbl, "tests_kruskal_imv_hours_by_exposure_tertile")

# 4) “Max support within 24h” vs exposure tertile (exact test on collapsed categories)
# Collapse to IMV vs not-IMV for stability
tab_max24_imv_pm25 <- dat2 %>%
  mutate(max24_imv = (as.character(max_support_24h) == "IMV")) %>%
  filter(!is.na(pm25_tertile), !is.na(max24_imv)) %>%
  with(table(pm25_tertile, max24_imv))

fisher_max24_pm25 <- broom::tidy(fisher.test(tab_max24_imv_pm25))
save_csv(as_tibble(tab_max24_imv_pm25, .name_repair="minimal"), "tab_max24_imv_by_pm25_tertile_raw")
save_csv(fisher_max24_pm25, "test_fisher_max24_imv_by_pm25_tertile")

# 5) Extubation <=48h among IMV-at-t0 patients vs exposure tertile (exact)
tab_extub_pm25 <- dat2 %>%
  filter(imv_t0) %>%
  mutate(extub48 = as.logical(extub_48h)) %>%
  filter(!is.na(pm25_tertile), !is.na(extub48)) %>%
  with(table(pm25_tertile, extub48))

if (all(dim(tab_extub_pm25) > 0)) {
  fisher_extub_pm25 <- broom::tidy(fisher.test(tab_extub_pm25))
  save_csv(as_tibble(tab_extub_pm25, .name_repair="minimal"), "tab_extub48_by_pm25_tertile_raw")
  save_csv(fisher_extub_pm25, "test_fisher_extub48_by_pm25_tertile")
} else {
  message("Skipping extub48 fisher test (empty table due to small n).")
}

# --------------------------- Plots (simple, interpretable) ---------------------------
# A) PM2.5 by tertile (jitter + box)
p_pm25_tertile <- dat2 %>%
  filter(!is.na(pm25_tertile)) %>%
  ggplot(aes(x=factor(pm25_tertile), y=pm25_1y)) +
  geom_boxplot() +
  geom_jitter(width=.12, alpha=.7) +
  labs(title="PM2.5 (1y pre-op) by tertile", x="PM2.5 tertile", y="PM2.5 (1y)") +
  theme_pub()
save_plot(p_pm25_tertile, "fig_pm25_by_tertile", w=7.5, h=5)

# B) IMV hours (total) by exposure tertile
p_imv_hours_pm25 <- dat2 %>%
  filter(!is.na(pm25_tertile)) %>%
  ggplot(aes(x=factor(pm25_tertile), y=imv_hours_total)) +
  geom_boxplot() +
  geom_jitter(width=.12, alpha=.7) +
  labs(title="IMV hours during ICU stay by PM2.5 tertile", x="PM2.5 tertile", y="IMV hours (ICU total)") +
  theme_pub()
save_plot(p_imv_hours_pm25, "fig_imv_hours_total_by_pm25_tertile", w=8.5, h=5)

# C) Max support within 24h distribution (stacked, as proportions)
rs_levels <- c("Room Air","Low-flow O2","HFNC","NIV/CPAP","IMV","None/Other")
p_max24 <- dat2 %>%
  mutate(max_support_24h = factor(as.character(max_support_24h), levels=rs_levels)) %>%
  filter(!is.na(pm25_tertile)) %>%
  count(pm25_tertile, max_support_24h) %>%
  group_by(pm25_tertile) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=factor(pm25_tertile), y=prop, fill=max_support_24h)) +
  geom_col(position="fill") +
  scale_y_continuous(labels=percent_format()) +
  labs(title="Max respiratory support within 24h by PM2.5 tertile",
       x="PM2.5 tertile", y="Proportion", fill="Max support (24h)") +
  theme_pub()
save_plot(p_max24, "fig_max_support_24h_by_pm25_tertile", w=9, h=5.5)

# D) Simple respiratory course label (for quick reviewer narrative)
dat2 <- dat2 %>%
  mutate(
    resp_course = case_when(
      imv_t0 & extub_24h ~ "IMV rapid extub ≤24h",
      imv_t0 & !extub_24h & extub_48h ~ "IMV extub 24–48h",
      imv_t0 & !extub_48h ~ "IMV prolonged/unclear",
      !imv_t0 & imv_any_72 ~ "Escalated to IMV ≤72h",
      !imv_t0 & !imv_any_72 ~ "No IMV in 72h",
      TRUE ~ "Other/unknown"
    )
  )

p_course_pm25 <- dat2 %>%
  ggplot(aes(x=resp_course, y=pm25_1y)) +
  geom_boxplot() +
  geom_jitter(width=.12, alpha=.6) +
  labs(title="PM2.5 (1y) by respiratory course", x=NULL, y="PM2.5 (1y)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle=25, hjust=1))
save_plot(p_course_pm25, "fig_pm25_by_resp_course", w=11, h=5.5)

# --------------------------- Federated summary package ---------------------------
# Core site-level summary (single row) for pooling across sites
site_summary <- dat2 %>%
  summarise(
    site_name = SITE_NAME,
    n = n(),
    n_death_in_hosp = sum(death_in_hosp, na.rm=TRUE),
    n_death_or_hospice = sum(death_or_hospice, na.rm=TRUE),
    n_imv_t0 = sum(imv_t0, na.rm=TRUE),
    n_imv_any_72 = sum(imv_any_72, na.rm=TRUE),
    median_icu_los_days = median(icu_los_days, na.rm=TRUE),
    median_imv_hours_total = median(imv_hours_total, na.rm=TRUE),
    median_imv_hours_72 = median(imv_hours_72, na.rm=TRUE),
    pm25_1y_median = median(pm25_1y, na.rm=TRUE),
    no2_1y_median  = median(no2_1y,  na.rm=TRUE)
  )

save_csv(site_summary, "federated_site_summary")

# Additionally save compact “2xK” count tables for meta-analysis
# (These avoid sharing individual-level data if needed.)
federated_counts <- list(
  death_by_pm25_tertile = as.data.frame.matrix(table(dat2$pm25_tertile, dat2$death_in_hosp)),
  death_by_no2_tertile  = as.data.frame.matrix(table(dat2$no2_tertile,  dat2$death_in_hosp)),
  max24_imv_by_pm25_tertile = as.data.frame.matrix(tab_max24_imv_pm25)
)

# Write each as CSV
for (nm in names(federated_counts)) {
  df <- federated_counts[[nm]] %>%
    tibble::rownames_to_column("exposure_level")
  save_csv(df, paste0("federated_counts_", nm))
}
















message("\nDONE. All outputs saved to:\n  ", out_dir, "\n")




