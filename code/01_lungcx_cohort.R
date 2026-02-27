# ================================================================================================
# CLIF | Lung Cancer POA ICU Cohort + Exposome + Respiratory Support Course Features (Federated)
# PI: Peter Graffy
#
# Pivot goal:
#   1) Build cohort: Lung cancer diagnosis present-on-arrival (POA=1) + any ICU segment
#   2) Anchor t0 at first ICU-in
#   3) Link county-year PM2.5 / NO2 exposures by (county_code, admission_year)
#   4) Derive respiratory support course features (ICU-wide + first 72h) using hour-binned aggregates
#   5) Save federated-friendly outputs (no raw timestamps)
#
# Notes:
#   - Assumes CLIF tables exist under tables_path and are named (or can be normalized to) clif_<table>
#   - Avoids procedure-code reliance for primary cohort build (resection can be added as optional tag later)
#   - Outputs to repo/output/run_[SITE]_[DATE]/
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

ADULT_AGE_YEARS <- 18
T0_ANCHOR       <- "icu_in"  # fixed for this analysis

# dx prefixes live here (ICD9/10-CM lung cancer prefixes)
CODES_PATH <- file.path(repo, "outlier-thresholds", "table.tsv")
stopifnot(file.exists(CODES_PATH))

# exposome input files (county-year)
EXPO_DIR  <- file.path(repo, "exposome")
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

required_raw <- c("patient","hospitalization","adt","hospital_diagnosis","respiratory_support")
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
    has_geo  = !is.na(census_tract) | !is.na(county_code)
  )

# B) Lung cancer dx POA=1
lung_poa <- hospital_dx %>%
  filter(poa_present == 1) %>%
  filter(code_matches_any_prefix(icd_code, lung_dx_prefixes)) %>%
  distinct(hospitalization_id) %>%
  mutate(has_lung_cancer_poa = TRUE)

# C) ICU segments
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

# Steps for flow
step0 <- base
step1 <- step0 %>% filter(adult)
step2 <- step1 %>% filter(has_demo)
step3 <- step2 %>% filter(has_geo)
step4 <- step3 %>% inner_join(lung_poa, by="hospitalization_id")
step5 <- step4 %>% inner_join(icu_segments %>% distinct(hospitalization_id), by="hospitalization_id")

# Index ICU = first ICU segment
index_icu <- icu_segments %>%
  group_by(hospitalization_id) %>%
  slice_min(icu_in_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(hospitalization_id, icu_in_time, icu_out_time)

cohort_lung <- step5 %>%
  inner_join(index_icu, by="hospitalization_id") %>%
  left_join(icu_bounds, by="hospitalization_id") %>%
  mutate(t0 = icu_in_time) %>%
  transmute(
    cohort = "lung_cancer_poa_icu",
    patient_id, hospitalization_id,
    admission_dttm, discharge_dttm,
    discharge_category, death_in_hosp, hospice_discharge, death_or_hospice,
    first_icu_in, last_icu_out, icu_los_hours,
    age_years, sex_category, race_category, ethnicity_category,
    census_tract, county_code, zipcode_five_digit, zipcode_nine_digit,
    icu_in_time, icu_out_time,
    t0,
    adult, has_demo, has_geo
  )

flow_lung <- tibble(
  step = c(
    "Hospitalizations in date window",
    glue(">= {ADULT_AGE_YEARS} years"),
    "Demographics present",
    "Geography present",
    "Lung cancer dx POA present",
    "Any ICU segment present"
  ),
  remaining = c(nrow(step0), nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5))
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
  left_join(icu_segments %>% distinct(hospitalization_id) %>% mutate(icu_ok=TRUE), by="hospitalization_id") %>%
  mutate(
    reason = coalesce(
      reason,
      ifelse(is.na(dx_ok), "No lung cancer dx POA=1 matching prefixes", NA_character_),
      ifelse(is.na(icu_ok), "No ICU stay (no ADT segment with ICU location_category)", NA_character_),
      "Other"
    )
  ) %>%
  filter(!(hospitalization_id %in% cohort_lung$hospitalization_id)) %>%
  select(patient_id, hospitalization_id, reason)

cat("\nCohort selection summary:\n")
cat("  Base hospitalizations:         ", nrow(step0), "\n", sep="")
cat("  Lung cancer POA eligible:      ", nrow(step4), "\n", sep="")
cat("  Eligible + ICU:                ", nrow(step5), "\n", sep="")
cat("  Final cohort included:         ", nrow(cohort_lung), "\n", sep="")


# --------------------------- Respiratory support course features ---------------------------
icu_windows <- cohort_lung %>%
  transmute(hospitalization_id,
            t0,
            icu_start = first_icu_in,
            icu_end   = last_icu_out)

rs_cohort <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(icu_windows, by = "hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    in_icu = recorded_dttm >= icu_start & recorded_dttm <= icu_end
  ) %>%
  filter(in_icu)

# Device intensity ranking (EDIT as needed for your site value set)
device_rank <- function(x) {
  x <- toupper(trimws(as.character(x)))
  case_when(
    x %in% c("ECMO") ~ 5,
    x %in% c("IMV","INVASIVE","MECHANICALVENTILATION") ~ 4,
    x %in% c("NIV","BIPAP","CPAP") ~ 3,
    x %in% c("HFNC","HIGHFLOW") ~ 2,
    x %in% c("NC","NASALCANNULA","O2","OXYGEN") ~ 1,
    is.na(x) | x == "" ~ NA_real_,
    TRUE ~ 0
  )
}

rs_cohort <- rs_cohort %>% mutate(rank = device_rank(device_category))

# A) First device at/after t0 (within +6h window)
first_device <- rs_cohort %>%
  filter(dt_h >= 0, dt_h <= 6) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(hospitalization_id, device_t0 = device_category)

# B) IMV hours in ICU and first 72h (hour-binned)
imv_hours <- rs_cohort %>%
  mutate(
    dt_h_floor = floor(dt_h),
    is_imv = toupper(device_category) %in% c("IMV","INVASIVE","MECHANICALVENTILATION")
  ) %>%
  filter(!is.na(dt_h_floor)) %>%
  group_by(hospitalization_id) %>%
  summarize(
    imv_hours_icu = n_distinct(dt_h_floor[is_imv]),
    imv_hours_72h = n_distinct(dt_h_floor[is_imv & dt_h_floor >= 0 & dt_h_floor <= 72]),
    any_imv_icu   = any(is_imv, na.rm = TRUE),
    any_imv_72h   = any(is_imv & dt_h_floor >= 0 & dt_h_floor <= 72, na.rm = TRUE),
    .groups="drop"
  )

# C) Escalation/de-escalation counts using max rank per hour
rs_hourly_state <- rs_cohort %>%
  filter(!is.na(rank), !is.na(dt_h)) %>%
  mutate(dt_h_floor = floor(dt_h)) %>%
  group_by(hospitalization_id, dt_h_floor) %>%
  summarize(
    rank_h = max(rank, na.rm = TRUE),
    .groups="drop"
  ) %>%
  arrange(hospitalization_id, dt_h_floor) %>%
  group_by(hospitalization_id) %>%
  mutate(
    delta = rank_h - lag(rank_h),
    escalation = ifelse(!is.na(delta) & delta > 0, 1L, 0L),
    deescalation = ifelse(!is.na(delta) & delta < 0, 1L, 0L)
  ) %>%
  summarize(
    escalations_icu = sum(escalation, na.rm = TRUE),
    deescalations_icu = sum(deescalation, na.rm = TRUE),
    .groups="drop"
  )

# D) Combine course features
rs_features <- cohort_lung %>%
  transmute(hospitalization_id, t0) %>%
  left_join(first_device, by="hospitalization_id") %>%
  left_join(imv_hours, by="hospitalization_id") %>%
  left_join(rs_hourly_state, by="hospitalization_id") %>%
  mutate(across(c(imv_hours_icu, imv_hours_72h, escalations_icu, deescalations_icu), ~coalesce(.x, 0))) %>%
  mutate(
    any_imv_icu = coalesce(any_imv_icu, FALSE),
    any_imv_72h = coalesce(any_imv_72h, FALSE)
  )

# --------------------------- Exposome linkage (county-year) ---------------------------
pm25 <- readr::read_csv(PM25_PATH, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  rename(county = geoid, pm25 = pm25_mean) %>%
  transmute(
    county = norm_county(county),
    year   = as.integer(year),
    pm25   = as.numeric(pm25)
  ) %>%
  distinct(county, year, .keep_all = TRUE)

no2 <- readr::read_csv(NO2_PATH, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  rename(county = geoid, no2 = no2_mean) %>%
  transmute(
    county = norm_county(county),
    year   = as.integer(year),
    no2    = as.numeric(no2)
  ) %>%
  distinct(county, year, .keep_all = TRUE)

expo <- pm25 %>%
  full_join(no2, by = c("county","year")) %>%
  arrange(county, year)

# Helper: trailing mean requiring full window coverage
trailing_mean_full <- function(x, k) {
  # x is ordered by year within county
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (i >= k) {
      w <- x[(i - k + 1):i]
      if (all(!is.na(w))) out[i] <- mean(w)
    }
  }
  out
}

expo_roll <- expo %>%
  group_by(county) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    pm25_1y = pm25,
    pm25_3y = trailing_mean_full(pm25, 3),
    pm25_5y = trailing_mean_full(pm25, 5),
    no2_1y  = no2,
    no2_3y  = trailing_mean_full(no2, 3),
    no2_5y  = trailing_mean_full(no2, 5)
  ) %>%
  ungroup() %>%
  select(county, year, pm25_1y, pm25_3y, pm25_5y, no2_1y, no2_3y, no2_5y)

# Join to cohort (and keep only the new exposure columns if you want)
cohort_expo <- cohort_lung %>%
  mutate(
    county = norm_county(county_code),
    year   = as.integer(lubridate::year(as.Date(admission_dttm)))
  ) %>%
  left_join(expo_roll, by = c("county","year")) %>%
  select(-county, -year)

# QC
cohort_expo %>%
  summarize(
    pm25_3y_unique = n_distinct(pm25_3y[!is.na(pm25_3y)]),
    pm25_5y_unique = n_distinct(pm25_5y[!is.na(pm25_5y)]),
    no2_3y_unique  = n_distinct(no2_3y[!is.na(no2_3y)]),
    no2_5y_unique  = n_distinct(no2_5y[!is.na(no2_5y)]),
    pm25_3y_missing = sum(is.na(pm25_3y)),
    pm25_5y_missing = sum(is.na(pm25_5y)),
    no2_3y_missing  = sum(is.na(no2_3y)),
    no2_5y_missing  = sum(is.na(no2_5y))
  ) %>% print()


# --------------------------- Final analysis-ready bundle (federated friendly) ---------------------------
analysis_ready <- cohort_lung %>%
  select(
    hospitalization_id,
    admission_dttm, discharge_dttm,
    death_in_hosp, hospice_discharge, death_or_hospice,
    icu_los_hours,
    age_years, sex_category, race_category, ethnicity_category,
    census_tract, county_code
  ) %>%
  left_join(rs_features %>% select(-t0), by = "hospitalization_id") %>%
  left_join(
    cohort_expo %>%
      select(hospitalization_id, pm25_1y, pm25_3y, pm25_5y, no2_1y, no2_3y, no2_5y),
    by = "hospitalization_id"
  )

cat("\nDone. Outputs in: ", out_dir, "\n", sep="")






