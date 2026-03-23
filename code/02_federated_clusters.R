# ================================================================================================
# CLIF | Lung Cancer ICU Trajectory Phenotypes x Air Pollution Association
#
# PURPOSE
#   This script:
#     1) Builds hourly (0–72h) respiratory support and ARF subtype panels
#     2) Derives hourly vasoactive use (for static summaries; not used in final clustering)
#     3) Collapses respiratory support into clinically interpretable tiers
#     4) Performs state-sequence clustering on: Collapsed RS | ARF subtype
#     5) Characterizes each cluster using:
#          - demographics / outcomes
#          - vasoactive use
#          - Charlson
#          - advanced cancer proxies
#          - first-24h SOFA
#     6) Runs air pollution association models:
#          - multinomial cluster assignment
#          - mortality
#          - interaction (NO2 × cluster, PM2.5 × cluster)
#          - ICU LOS
#          - ordinal trajectory severity
#          - 72h landmark Cox survival
#          - mediation-style attenuation analyses
#     7) Exports publication-ready figures and standardized model output tables
#     8) Exports federated-friendly cluster centroids / static summaries / trajectory summaries
#
# ASSUMPTIONS
#   The following objects already exist in memory from your upstream cohort builder:
#     - cohort_lung
#     - analysis_ready   (with outcomes + exposures already joined: pm25_1y/3y/5y, no2_1y/3y/5y)
#     - rs_raw           (CLIF respiratory support)
#     - icu_segments
#     - hospital_dx
#     - clif_tables      (named list of CLIF tables)
#     - safe_posix(), as_num(), norm_code(), code_matches_any_prefix()
#     - get_min(), save_csv(), save_plot(), make_name(), out_dir, theme_pub()
#     - device_state_clif(), device_rank_clif()  # exact CLIF mappings defined upstream
#
# NOTES
#   - Final clustering uses Collapsed RS | ARF subtype (not vaso), per latest analysis decisions
#   - Vasoactive use is retained for cluster characterization / static summaries
#   - All output file names are standardized for site-level federated export
# ================================================================================================

# Source 01 to build cohort (needed when running via Rscript in a fresh process)
if (!exists("cohort_lung")) {
  source("code/01_lungcx_cohort.R")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(TraMineR)
  library(cluster)
  library(data.table)
  library(nnet)
  library(MASS)
  library(survival)
  library(comorbidity)
  library(ggplot2)
  library(sf)
  library(tigris)
})

# ---------------------------
# Global params
# ---------------------------
H_MAX <- 72L
ROOM_AIR_FIO2 <- 0.21
JOIN_NEAR_H   <- 1
HYPERPAIR_H   <- 2
CARRY_H       <- 6L
K_RA          <- 5L   # set from silhouette / interpretability review

# ---------------------------
# Guardrails
# ---------------------------
stopifnot("hospitalization_id" %in% names(cohort_lung))
stopifnot("t0" %in% names(cohort_lung))
stopifnot("pm25_5y" %in% names(analysis_ready))
stopifnot("no2_5y" %in% names(analysis_ready))
stopifnot("death_or_hospice" %in% names(analysis_ready))
stopifnot("icu_los_hours" %in% names(analysis_ready))

# ---------------------------
# Small helpers
# ---------------------------
wmean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

carry_forward <- function(x, carry = 6L) {
  out <- rep(FALSE, length(x))
  idx <- which(x)
  for (i in idx) out[i:min(length(x), i + carry)] <- TRUE
  out
}

collapse_rs <- function(rs_state) {
  x <- toupper(trimws(as.character(rs_state)))
  case_when(
    x %in% c("ROOM AIR") ~ "ROOM AIR",
    x %in% c("NASAL CANNULA", "FACE MASK", "TRACH COLLAR") ~ "LOW_O2",
    x %in% c("HIGH FLOW NC", "CPAP", "NIPPV") ~ "NIV",
    x %in% c("IMV") ~ "IMV",
    x %in% c("OTHER", "UNK", "", NA_character_) ~ "OTHER",
    TRUE ~ "OTHER"
  )
}

# standardized model-saving helpers
save_model_table <- function(df, stub) {
  df <- df %>% mutate(across(everything(), ~ ifelse(is.infinite(.x), NA, .x)))
  save_csv(df, stub)
}

tidy_glm_model <- function(model, model_name) {
  s <- summary(model)
  cf <- as.data.frame(s$coefficients)
  term <- rownames(cf)
  out <- tibble(
    model = model_name,
    term = term,
    estimate = cf[[1]],
    std_error = cf[[2]],
    statistic = cf[[3]],
    p_value = cf[[4]]
  )
  if (inherits(model, "glm") && family(model)$family == "binomial") {
    out <- out %>%
      mutate(
        effect = exp(estimate),
        conf_low = exp(estimate - 1.96 * std_error),
        conf_high = exp(estimate + 1.96 * std_error),
        effect_type = "OR"
      )
  } else if (inherits(model, "glm") && grepl("poisson", family(model)$family)) {
    out <- out %>%
      mutate(
        effect = exp(estimate),
        conf_low = exp(estimate - 1.96 * std_error),
        conf_high = exp(estimate + 1.96 * std_error),
        effect_type = "RR_like"
      )
  } else {
    out <- out %>%
      mutate(effect = NA_real_, conf_low = NA_real_, conf_high = NA_real_, effect_type = NA_character_)
  }
  out
}

tidy_multinom_model <- function(model, model_name) {
  s <- summary(model)
  coefs <- s$coefficients
  ses <- s$standard.errors
  rn <- rownames(coefs)
  cn <- colnames(coefs)
  expand_grid(outcome_level = rn, term = cn) %>%
    mutate(
      estimate = map2_dbl(outcome_level, term, ~ coefs[.x, .y]),
      std_error = map2_dbl(outcome_level, term, ~ ses[.x, .y]),
      statistic = estimate / std_error,
      p_value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
      effect = exp(estimate),
      conf_low = exp(estimate - 1.96 * std_error),
      conf_high = exp(estimate + 1.96 * std_error),
      effect_type = "RRR",
      model = model_name
    ) %>%
    relocate(model, outcome_level, term)
}

tidy_polr_model <- function(model, model_name) {
  cf <- coef(summary(model))
  tibble(
    model = model_name,
    term = rownames(cf),
    estimate = cf[, "Value"],
    std_error = cf[, "Std. Error"],
    statistic = cf[, "t value"],
    p_value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
    effect = exp(estimate),
    conf_low = exp(estimate - 1.96 * std_error),
    conf_high = exp(estimate + 1.96 * std_error),
    effect_type = "Ordinal_OR"
  )
}

tidy_cox_model <- function(model, model_name) {
  s <- summary(model)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  tibble(
    model = model_name,
    term = rownames(cf),
    estimate = cf[["coef"]],
    std_error = cf[["se(coef)"]],
    statistic = cf[["z"]],
    p_value = cf[["Pr(>|z|)"]],
    effect = ci[["exp(coef)"]],
    conf_low = ci[["lower .95"]],
    conf_high = ci[["upper .95"]],
    effect_type = "HR"
  )
}

pick_k_silhouette <- function(dist_mat, k_min = 2, k_max = 10, hclust_method = "ward.D2") {
  d <- if (inherits(dist_mat, "dist")) dist_mat else as.dist(dist_mat)
  hc <- hclust(d, method = hclust_method)
  res <- map_dfr(k_min:k_max, function(k) {
    cl <- cutree(hc, k = k)
    sil <- silhouette(cl, d)
    tibble(
      k = k,
      avg_sil_width = mean(sil[, "sil_width"], na.rm = TRUE),
      min_sil_width = min(sil[, "sil_width"], na.rm = TRUE)
    )
  })
  p <- ggplot(res, aes(k, avg_sil_width)) +
    geom_line() +
    geom_point() +
    labs(title = "Average silhouette width vs K", x = "K", y = "Average silhouette width") +
    theme_pub()
  list(table = res, plot = p, hc = hc)
}

state_colors <- function(state) {
  parts <- strsplit(state, "\\|")[[1]]
  rs  <- parts[1]
  arf <- parts[2]
  
  base <- case_when(
    rs == "ROOM AIR" ~ "#2ca25f",
    rs == "LOW_O2"   ~ "#e6ab02",
    rs == "NIV"      ~ "#3182bd",
    rs == "IMV"      ~ "#de2d26",
    rs == "OTHER"    ~ "#636363",
    TRUE ~ "#969696"
  )
  
  adjust <- function(col, fac) {
    rgb_col <- col2rgb(col) / 255
    rgb_adj <- pmin(1, rgb_col * fac)
    rgb(rgb_adj[1], rgb_adj[2], rgb_adj[3])
  }
  
  case_when(
    arf == "NO_ARF"    ~ adjust(base, 1.2),
    arf == "ARF_HYPOX" ~ adjust(base, 1.0),
    arf == "ARF_HYPER" ~ adjust(base, 0.8),
    arf == "ARF_MIXED" ~ adjust(base, 0.6),
    TRUE ~ base
  )
}

# ---------------------------
# Keys / windows
# ---------------------------
keys <- cohort_lung %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    t0 = as.POSIXct(t0, tz = "UTC")
  )

win72 <- keys %>%
  transmute(
    hospitalization_id,
    win_start = t0,
    win_end   = t0 + dhours(H_MAX)
  )

# ================================================================================================
# 1) Hourly Respiratory Support Panel
# ================================================================================================

rs_for_traj <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(keys, by = "hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units = "hours")),
    h = floor(dt_h),
    rs_state = device_state_clif(device_category),
    rank = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

hourly_rs <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
 dplyr::select(hospitalization_id, h, rs_state)

rs_panel <- hourly_rs %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(rs_state, .direction = "down") %>%
  tidyr::fill(rs_state, .direction = "up") %>%
  mutate(
    rs_state = coalesce(rs_state, "UNK"),
    rs_state_collapsed = collapse_rs(rs_state)
  ) %>%
  ungroup()

#save_csv(rs_panel, "panel_rs_hourly_collapsed")

# ================================================================================================
# 2) Hourly Vasoactive Use (static summaries only; not used in final clustering)
# ================================================================================================

med_cont <- get_min("medication_admin_continuous", c(
  "hospitalization_id", "admin_dttm", "med_group", "mar_action_group"
)) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admin_dttm = safe_posix(admin_dttm),
    med_group = tolower(as.character(med_group)),
    mar_action_group = tolower(as.character(mar_action_group))
  )

vaso_hourly <- med_cont %>%
  filter(
    med_group == "vasoactives",
    mar_action_group == "administered",
    !is.na(admin_dttm)
  ) %>%
  inner_join(keys, by = "hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(admin_dttm, t0, units = "hours")),
    h = floor(dt_h)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX) %>%
  distinct(hospitalization_id, h) %>%
  mutate(vaso_h = 1L)

vaso_static <- tibble(hospitalization_id = unique(keys$hospitalization_id)) %>%
  left_join(
    vaso_hourly %>%
      group_by(hospitalization_id) %>%
      summarize(
        vaso_any_72h = as.integer(any(vaso_h == 1L, na.rm = TRUE)),
        vaso_hours_72h = sum(vaso_h, na.rm = TRUE),
        vaso_mean_72h = mean(vaso_h, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "hospitalization_id"
  ) %>%
  mutate(
    vaso_any_72h = coalesce(vaso_any_72h, 0L),
    vaso_hours_72h = coalesce(vaso_hours_72h, 0),
    vaso_mean_72h = coalesce(vaso_mean_72h, 0)
  )

#save_csv(vaso_static, "vaso_static_72h")

# ================================================================================================
# 3) Hourly ARF subtype panel
# ================================================================================================

vitals <- clif_tables[["clif_vitals"]]
labs   <- clif_tables[["clif_labs"]]

spo2_win <- vitals %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    vital_category = tolower(as.character(vital_category)),
    vital_value = as_num(vital_value)
  ) %>%
  filter(vital_category == "spo2") %>%
  inner_join(win72, by = "hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, spo2_time = recorded_dttm, spo2 = vital_value)

labs_win <- labs %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm = safe_posix(lab_result_dttm),
    lab_category = tolower(as.character(lab_category)),
    lab_value_numeric = as_num(lab_value_numeric)
  ) %>%
  filter(lab_category %in% c("po2_arterial", "pco2_arterial", "ph_arterial")) %>%
  inner_join(win72, by = "hospitalization_id") %>%
  filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
  transmute(hospitalization_id, lab_time = lab_result_dttm, lab_category, val = lab_value_numeric)

fio2_win <- rs_raw %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    fio2_set = as_num(fio2_set)
  ) %>%
  inner_join(win72, by = "hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, fio2_time = recorded_dttm, fio2_set)

setDT(spo2_win); setDT(fio2_win); setDT(labs_win)
setkey(spo2_win, hospitalization_id, spo2_time)
setkey(fio2_win, hospitalization_id, fio2_time)

spo2_win[, spo2_time_keep := spo2_time]
spo2_fio2 <- fio2_win[
  spo2_win, roll = "nearest", on = .(hospitalization_id, fio2_time = spo2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(spo2_time_keep, fio2_time, units = "hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, t = spo2_time_keep, spo2, fio2_set,
    on_room_air = !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6)
]

po2dt <- labs_win[lab_category == "po2_arterial", .(hospitalization_id, po2_time = lab_time, po2 = val)]
setkey(po2dt, hospitalization_id, po2_time)
po2dt[, po2_time_keep := po2_time]
setkey(fio2_win, hospitalization_id, fio2_time)

po2_fio2 <- fio2_win[
  po2dt, roll = "nearest", on = .(hospitalization_id, fio2_time = po2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(po2_time_keep, fio2_time, units = "hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, t = po2_time_keep, po2, fio2_set,
    pf_ratio = fifelse(!is.na(fio2_set) & fio2_set > 0, po2 / fio2_set, as.numeric(NA)))
]

pco2dt <- labs_win[lab_category == "pco2_arterial", .(hospitalization_id, pco2_time = lab_time, pco2 = val)]
phdt   <- labs_win[lab_category == "ph_arterial", .(hospitalization_id, ph_time = lab_time, ph = val)]
setkey(pco2dt, hospitalization_id, pco2_time)
setkey(phdt, hospitalization_id, ph_time)

pco2dt[, pco2_time_keep := pco2_time]
phdt[, ph_time_keep := ph_time]

hyper_pairs <- phdt[
  pco2dt, roll = "nearest", on = .(hospitalization_id, ph_time = pco2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(pco2_time_keep, ph_time, units = "hours")))
][
  timediff_h <= as.numeric(HYPERPAIR_H),
  .(hospitalization_id, t = pco2_time_keep, pco2, ph,
    hyper_hit = (pco2 >= 45 & ph < 7.35))
]

keys_dt <- as.data.table(keys)
setkey(keys_dt, hospitalization_id)

to_hour_bin <- function(dt, time_col) {
  dt[keys_dt, on = "hospitalization_id"][,
                                         h := floor(as.numeric(difftime(get(time_col), t0, units = "hours")))
  ][h >= 0 & h <= H_MAX]
}

spo2_fio2_h <- to_hour_bin(spo2_fio2, "t")[,
                                           .(hypox_roomair_spo2 = any(spo2 < 90 & on_room_air, na.rm = TRUE)),
                                           by = .(hospitalization_id, h)
]

po2_fio2_h <- to_hour_bin(po2_fio2, "t")[,
                                         .(
                                           hypox_roomair_po2 = any(po2 <= 60 & !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6, na.rm = TRUE),
                                           hypox_pf = any(!is.na(pf_ratio) & pf_ratio <= 300, na.rm = TRUE)
                                         ),
                                         by = .(hospitalization_id, h)
]

hyper_h <- to_hour_bin(hyper_pairs, "t")[,
                                         .(hypercapnic = any(hyper_hit, na.rm = TRUE)),
                                         by = .(hospitalization_id, h)
]

arf_hourly <- Reduce(
  function(x, y) merge(x, y, by = c("hospitalization_id", "h"), all = TRUE),
  list(spo2_fio2_h, po2_fio2_h, hyper_h)
) %>%
  as_tibble() %>%
  mutate(
    hypox_roomair_spo2 = coalesce(hypox_roomair_spo2, FALSE),
    hypox_roomair_po2  = coalesce(hypox_roomair_po2, FALSE),
    hypox_pf           = coalesce(hypox_pf, FALSE),
    hypercapnic        = coalesce(hypercapnic, FALSE),
    hypox_h = hypox_roomair_spo2 | hypox_roomair_po2 | hypox_pf,
    hyper_h = hypercapnic
  ) %>%
  dplyr::select(hospitalization_id, h, hypox_h, hyper_h) %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX, fill = list(hypox_h = FALSE, hyper_h = FALSE)) %>%
  arrange(h) %>%
  mutate(
    hypox_h = carry_forward(hypox_h, carry = CARRY_H),
    hyper_h = carry_forward(hyper_h, carry = CARRY_H),
    arf_subtype_h = case_when(
      hypox_h & hyper_h ~ "ARF_MIXED",
      hypox_h & !hyper_h ~ "ARF_HYPOX",
      !hypox_h & hyper_h ~ "ARF_HYPER",
      TRUE ~ "NO_ARF"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(hospitalization_id, h, arf_subtype_h)

#save_csv(arf_hourly, "panel_arf_hourly")

# ================================================================================================
# 4) Final clustering panel: Collapsed RS | ARF subtype
# ================================================================================================

panel_ra <- rs_panel %>%
  left_join(arf_hourly, by = c("hospitalization_id", "h")) %>%
  mutate(
    arf_subtype_h = coalesce(arf_subtype_h, "NO_ARF"),
    rs_state_collapsed = coalesce(rs_state_collapsed, "OTHER"),
    state_ra = paste(rs_state_collapsed, arf_subtype_h, sep = "|")
  )

#save_csv(panel_ra, "panel_ra_hourly")

# ================================================================================================
# 5) Sequence clustering: Collapsed RS | ARF subtype
# ================================================================================================

seq_wide_ra <- panel_ra %>%
  mutate(h_col = paste0("H", h)) %>%
  dplyr::select(hospitalization_id, h_col, state_ra) %>%
  pivot_wider(names_from = h_col, values_from = state_ra)

seq_cols_ra <- seq_wide_ra %>% dplyr::select(matches("^H\\d+$"))
alphabet_ra <- sort(unique(unlist(seq_cols_ra)))
cpal_ra <- sapply(alphabet_ra, state_colors)

seq_obj_ra <- TraMineR::seqdef(
  seq_cols_ra,
  alphabet = alphabet_ra,
  states   = alphabet_ra,
  labels   = alphabet_ra,
  cpal     = cpal_ra,
  xtstep   = 6
)

dist_ra <- TraMineR::seqdist(seq_obj_ra, method = "OM", indel = 1, sm = "TRATE")
hc_ra <- hclust(as.dist(dist_ra), method = "ward.D2")

sil_out <- pick_k_silhouette(dist_ra, k_min = 2, k_max = 10)
save_csv(sil_out$table, "cluster_silhouette_ra")
save_plot(sil_out$plot, "cluster_silhouette_ra", w = 8, h = 5)

cluster_ra <- cutree(hc_ra, k = K_RA)

traj_assign_ra <- tibble(
  hospitalization_id = seq_wide_ra$hospitalization_id,
  traj_cluster_ra = factor(as.integer(cluster_ra))
)


png(file.path(out_dir, make_name("traj_ra_seqrplot_discrete", "png")), width = 1800, height = 2600, res = 160)
TraMineR::seqrplot(
  seq_obj_ra,
  group = cluster_ra,
  diss  = dist_ra,
  main = "Representative sequences (Collapsed RS | ARF subtype)"
)
dev.off()

png(file.path(out_dir, make_name("traj_ra_seqdplot_by_cluster", "png")), width = 1800, height = 2200, res = 160)
TraMineR::seqdplot(
  seq_obj_ra,
  group = cluster_ra,
  border = NA,
  main = "State distribution by cluster (Collapsed RS | ARF subtype)"
)
dev.off()

# ================================================================================================
# 6) Join trajectories to main analytic file
# ================================================================================================

analysis_ready <- analysis_ready %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  dplyr::select(-any_of("traj_cluster_ra")) %>%
  left_join(traj_assign_ra, by = "hospitalization_id")

# z-score exposures once
analysis_ready <- analysis_ready %>%
  mutate(
    pm25_5y_z = as.numeric(scale(pm25_5y)),
    no2_5y_z  = as.numeric(scale(no2_5y)),
    admit_year = year(admission_dttm),
    traj_cluster_ra = factor(traj_cluster_ra)
  )

# ================================================================================================
# 7) Charlson score
# ================================================================================================

dx_charlson <- hospital_dx %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    diagnosis_code = as.character(diagnosis_code),
    diagnosis_code_format = toupper(coalesce(as.character(diagnosis_code_format), "")),
    poa_present = suppressWarnings(as.integer(poa_present))
  ) %>%
  filter(
    poa_present == 1,
    !is.na(diagnosis_code),
    diagnosis_code != ""
  ) %>%
  mutate(
    code_clean = gsub("[^A-Za-z0-9]", "", diagnosis_code),
    icd_version = case_when(
      str_detect(diagnosis_code_format, "10") ~ "icd10",
      str_detect(diagnosis_code_format, "9")  ~ "icd9",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(icd_version))

charlson_icd10 <- dx_charlson %>%
  filter(icd_version == "icd10") %>%
  dplyr::select(id = hospitalization_id, code = code_clean) %>%
  comorbidity(id = "id", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)

charlson_icd10$charlson_score <- score(charlson_icd10, weights = "charlson", assign0 = FALSE)

charlson_all <- charlson_icd10 %>%
  dplyr::select(id, charlson_score) %>%
  group_by(id) %>%
  summarize(charlson_score = max(charlson_score, na.rm = TRUE), .groups = "drop") %>%
  rename(hospitalization_id = id)

analysis_ready <- analysis_ready %>%
  left_join(charlson_all, by = "hospitalization_id")

# ================================================================================================
# 8) Advanced cancer proxies
# ================================================================================================

metastatic_prefixes_icd10 <- c("C77", "C78", "C79", "C80")
metastatic_prefixes_icd9  <- c("196", "197", "198", "199")
pleural_effusion_icd10 <- c("J91", "C785")
pleural_effusion_icd9  <- c("51181")
cachexia_icd10 <- c("R64")
cachexia_icd9  <- c("7994")

advanced_cancer_flags <- hospital_dx %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    diagnosis_code = norm_code(diagnosis_code),
    diagnosis_code_format = toupper(coalesce(as.character(diagnosis_code_format), "")),
    poa_present = suppressWarnings(as.integer(poa_present)),
    icd10 = str_detect(diagnosis_code_format, "10"),
    icd9  = str_detect(diagnosis_code_format, "9")
  ) %>%
  filter(poa_present == 1) %>%
  group_by(hospitalization_id) %>%
  summarize(
    metastatic_cancer_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, metastatic_prefixes_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, metastatic_prefixes_icd9)),
      na.rm = TRUE
    ),
    malignant_pleural_effusion_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, pleural_effusion_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, pleural_effusion_icd9)),
      na.rm = TRUE
    ),
    cachexia_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, cachexia_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, cachexia_icd9)),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    advanced_cancer_any_poa =
      metastatic_cancer_poa |
      malignant_pleural_effusion_poa |
      cachexia_poa
  )

analysis_ready <- analysis_ready %>%
  left_join(advanced_cancer_flags, by = "hospitalization_id")

# ================================================================================================
# 9) First-24h SOFA
# ================================================================================================

source(file.path("utils", "sofa_calculator.R"))

icu_admit_times <- icu_segments %>%
  semi_join(cohort_lung, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>%
  summarize(
    icu_admit_time = min(icu_in_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    icu_admit_time = safe_posix(icu_admit_time)
  )

sofa_cohort <- cohort_lung %>%
  transmute(hospitalization_id = as.character(hospitalization_id)) %>%
  inner_join(icu_admit_times, by = "hospitalization_id")

vitals_df <- clif_tables[["clif_vitals"]] %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm)
  )

labs_df <- clif_tables[["clif_labs"]] %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm = safe_posix(lab_result_dttm)
  )

support_df <- rs_raw %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm)
  )

med_admin_df <- get_min("medication_admin_continuous", c(
  "hospitalization_id", "admin_dttm", "med_group", "med_category",
  "mar_action_group", "med_name", "med_dose", "med_dose_unit"
)) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admin_dttm = safe_posix(admin_dttm)
  )

scores_df <- clif_tables[["clif_patient_assessments"]] %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    assessment_category = as.character(assessment_category),
    numerical_value = as.numeric(numerical_value)
  ) %>%
  filter(!is.na(recorded_dttm), assessment_category == "gcs_total") %>%
  dplyr::select(hospitalization_id, recorded_dttm, assessment_category, numerical_value)

# safe_ts already defined in 01_lungcx_cohort.R

sofa_scores <- calculate_sofa(
  cohort_data = sofa_cohort,
  vitals_df = vitals_df,
  labs_df = labs_df,
  support_df = support_df,
  med_admin_df = med_admin_df,
  scores_df = scores_df,
  window_hours = 24,
  safe_ts = safe_ts
)

analysis_ready <- analysis_ready %>%
  dplyr::select(-any_of(c("sofa_total", "sofa_cv", "sofa_coag", "sofa_liver",
                           "sofa_renal", "sofa_resp", "sofa_cns"))) %>%
  left_join(
    sofa_scores %>%
      transmute(
        hospitalization_id = as.character(hospitalization_id),
        sofa_total, sofa_cv, sofa_coag, sofa_liver,
        sofa_renal, sofa_resp, sofa_cns
      ),
    by = "hospitalization_id"
  ) %>%
  mutate(across(starts_with("sofa_"), ~ coalesce(.x, 0)))

# ================================================================================================
# 10) Cluster descriptive tables, centroids, and trajectory summaries (federated outputs)
# ================================================================================================

# patient-level first-hit timing summaries
sig_panel <- rs_panel %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    h = as.integer(h),
    rs = as.character(rs_state_collapsed)
  ) %>%
  left_join(
    arf_hourly %>%
      transmute(
        hospitalization_id = as.character(hospitalization_id),
        h = as.integer(h),
        arf = as.character(arf_subtype_h)
      ),
    by = c("hospitalization_id", "h")
  ) %>%
  left_join(
    analysis_ready %>% transmute(hospitalization_id = as.character(hospitalization_id), traj_cluster_ra),
    by = "hospitalization_id"
  ) %>%
  mutate(
    rs  = toupper(coalesce(rs, "OTHER")),
    arf = toupper(coalesce(arf, "NO_ARF")),
    traj_cluster_ra = forcats::fct_drop(traj_cluster_ra)
  ) %>%
  filter(!is.na(traj_cluster_ra), h >= 0, h <= H_MAX)

first_hits <- sig_panel %>%
  group_by(hospitalization_id, traj_cluster_ra) %>%
  summarize(
    t_imv = suppressWarnings(min(h[rs == "IMV"], na.rm = TRUE)),
    t_arf = suppressWarnings(min(h[arf != "NO_ARF"], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    t_imv = ifelse(is.infinite(t_imv), NA, t_imv),
    t_arf = ifelse(is.infinite(t_arf), NA, t_arf)
  )

# static patient-level summaries for cluster characterization
cluster_static <- analysis_ready %>%
  left_join(vaso_static, by = "hospitalization_id") %>%
  left_join(first_hits %>% dplyr::select(hospitalization_id, t_imv, t_arf), by = "hospitalization_id") %>%
  filter(!is.na(traj_cluster_ra)) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    n = n(),
    age_mean = mean(age_years, na.rm = TRUE),
    age_sd = sd(age_years, na.rm = TRUE),
    male_pct = mean(sex_category == "Male", na.rm = TRUE),
    death_pct = mean(death_or_hospice, na.rm = TRUE),
    icu_los_median = median(icu_los_hours, na.rm = TRUE),
    imv72_pct = mean(any_imv_72h, na.rm = TRUE),
    imv72_median = median(imv_hours_72h, na.rm = TRUE),
    vaso_any_72h_pct = mean(vaso_any_72h, na.rm = TRUE),
    vaso_hours_72h_mean = mean(vaso_hours_72h, na.rm = TRUE),
    vaso_mean_72h = mean(vaso_mean_72h, na.rm = TRUE),
    pm25_5y_mean = mean(pm25_5y, na.rm = TRUE),
    no2_5y_mean = mean(no2_5y, na.rm = TRUE),
    charlson_mean = mean(charlson_score, na.rm = TRUE),
    charlson_median = median(charlson_score, na.rm = TRUE),
    sofa_total_mean = mean(sofa_total, na.rm = TRUE),
    sofa_resp_mean = mean(sofa_resp, na.rm = TRUE),
    sofa_cv_mean = mean(sofa_cv, na.rm = TRUE),
    sofa_renal_mean = mean(sofa_renal, na.rm = TRUE),
    sofa_liver_mean = mean(sofa_liver, na.rm = TRUE),
    sofa_coag_mean = mean(sofa_coag, na.rm = TRUE),
    sofa_cns_mean = mean(sofa_cns, na.rm = TRUE),
    metastatic_pct = mean(metastatic_cancer_poa, na.rm = TRUE),
    pleural_effusion_pct = mean(malignant_pleural_effusion_poa, na.rm = TRUE),
    cachexia_pct = mean(cachexia_poa, na.rm = TRUE),
    advanced_cancer_pct = mean(advanced_cancer_any_poa, na.rm = TRUE),
    t_imv_median = median(t_imv, na.rm = TRUE),
    t_arf_median = median(t_arf, na.rm = TRUE),
    .groups = "drop"
  )

save_csv(cluster_static, "cluster_static_summary")

# federated cluster centroids (compressed)
cluster_centroids <- sig_panel %>%
  left_join(vaso_hourly, by = c("hospitalization_id", "h")) %>%
  mutate(vaso_h = coalesce(vaso_h, 0L)) %>%
  group_by(traj_cluster_ra, hospitalization_id) %>%
  summarize(
    pct_room_air = mean(rs == "ROOM AIR", na.rm = TRUE),
    pct_low_o2   = mean(rs == "LOW_O2", na.rm = TRUE),
    pct_niv      = mean(rs == "NIV", na.rm = TRUE),
    pct_imv      = mean(rs == "IMV", na.rm = TRUE),
    pct_other    = mean(rs == "OTHER", na.rm = TRUE),
    pct_no_arf   = mean(arf == "NO_ARF", na.rm = TRUE),
    pct_hypox    = mean(arf == "ARF_HYPOX", na.rm = TRUE),
    pct_hyper    = mean(arf == "ARF_HYPER", na.rm = TRUE),
    pct_mixed    = mean(arf == "ARF_MIXED", na.rm = TRUE),
    pct_vaso     = mean(vaso_h == 1L, na.rm = TRUE),
    transition_count = sum(rs != lag(rs), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    across(-hospitalization_id, ~ mean(.x, na.rm = TRUE)),
    n_cluster = n(),
    .groups = "drop"
  )

save_csv(cluster_centroids, "cluster_centroids_federated")

# hourly trajectory summaries for federated harmonization
cluster_signature_rs <- sig_panel %>%
  count(traj_cluster_ra, h, rs) %>%
  group_by(traj_cluster_ra, h) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cluster_signature_arf <- sig_panel %>%
  count(traj_cluster_ra, h, arf) %>%
  group_by(traj_cluster_ra, h) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

save_csv(cluster_signature_rs, "cluster_signature_rs_hourly")
save_csv(cluster_signature_arf, "cluster_signature_arf_hourly")

# ================================================================================================
# 11) Cluster characterization figures
# ================================================================================================

p_sofa_total_box <- ggplot(analysis_ready, aes(traj_cluster_ra, sofa_total)) +
  geom_boxplot() +
  labs(
    x = "Trajectory cluster",
    y = "First 24h SOFA total",
    title = "Baseline ICU severity across respiratory trajectory clusters"
  ) +
  theme_pub()
save_plot(p_sofa_total_box, "cluster_sofa_total_boxplot", w = 8, h = 6)

p_sofa_cns_box <- ggplot(analysis_ready, aes(traj_cluster_ra, sofa_cns)) +
  geom_boxplot() +
  labs(
    x = "Trajectory cluster",
    y = "SOFA CNS score",
    title = "CNS organ dysfunction by trajectory cluster"
  ) +
  theme_pub()
save_plot(p_sofa_cns_box, "cluster_sofa_cns_boxplot", w = 8, h = 6)

cluster_domains <- analysis_ready %>%
  filter(!is.na(traj_cluster_ra)) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    resp = mean(sofa_resp, na.rm = TRUE),
    cv = mean(sofa_cv, na.rm = TRUE),
    renal = mean(sofa_renal, na.rm = TRUE),
    liver = mean(sofa_liver, na.rm = TRUE),
    coag = mean(sofa_coag, na.rm = TRUE),
    cns = mean(sofa_cns, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-traj_cluster_ra, names_to = "sofa_domain", values_to = "mean_score")

p_domains <- ggplot(cluster_domains, aes(sofa_domain, mean_score, fill = traj_cluster_ra)) +
  geom_col(position = "dodge") +
  labs(
    x = "SOFA domain",
    y = "Mean score",
    title = "Organ dysfunction signature by trajectory cluster"
  ) +
  theme_pub()
save_plot(p_domains, "cluster_sofa_domain_signature", w = 10, h = 6)

# trajectory signature plots
sig_rs <- cluster_signature_rs %>%
  mutate(rs = factor(rs, levels = c("ROOM AIR", "LOW_O2", "NIV", "IMV", "OTHER")))

p_rs <- ggplot(sig_rs, aes(x = h, y = prop, fill = rs)) +
  geom_area() +
  facet_wrap(~ traj_cluster_ra, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trajectory signatures by cluster: Respiratory support",
    x = "Hour from ICU t0",
    y = "Within-cluster proportion"
  ) +
  theme_pub()
save_plot(p_rs, "sig_rs_area_by_cluster", w = 12, h = 8)

sig_arf <- cluster_signature_arf %>%
  mutate(arf = factor(arf, levels = c("NO_ARF", "ARF_HYPOX", "ARF_HYPER", "ARF_MIXED")))

p_arf <- ggplot(sig_arf, aes(x = h, y = prop, fill = arf)) +
  geom_area() +
  facet_wrap(~ traj_cluster_ra, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trajectory signatures by cluster: ARF subtype",
    x = "Hour from ICU t0",
    y = "Within-cluster proportion"
  ) +
  theme_pub()
save_plot(p_arf, "sig_arf_area_by_cluster", w = 12, h = 8)

cum_curve <- function(df, time_var, label) {
  df %>%
    filter(!is.na(.data[[time_var]])) %>%
    count(traj_cluster_ra, t = .data[[time_var]]) %>%
    group_by(traj_cluster_ra) %>%
    arrange(t) %>%
    mutate(cum = cumsum(n) / sum(n), outcome = label) %>%
    ungroup()
}

cc_imv <- cum_curve(first_hits, "t_imv", "IMV start")
cc_arf <- cum_curve(first_hits, "t_arf", "Any ARF onset")
cc <- bind_rows(cc_imv, cc_arf)

p_cc <- ggplot(cc, aes(x = t, y = cum, color = traj_cluster_ra)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ outcome, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Timing signatures by cluster",
    x = "Hour from ICU t0",
    y = "Cumulative proportion with event",
    color = "Cluster"
  ) +
  theme_pub()
save_plot(p_cc, "sig_timing_cuminc_by_cluster", w = 10, h = 8)

# ================================================================================================
# 12) Cluster severity ranking
# ================================================================================================

cluster_severity <- analysis_ready %>%
  filter(!is.na(traj_cluster_ra)) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    mortality = mean(death_or_hospice, na.rm = TRUE),
    sofa = mean(sofa_total, na.rm = TRUE),
    imv_rate = mean(any_imv_72h, na.rm = TRUE),
    icu_los = median(icu_los_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    severity_rank = rank(mortality + sofa / 10 + imv_rate, ties.method = "first"),
    severity_rank = factor(severity_rank, ordered = TRUE)
  )

save_csv(cluster_severity, "cluster_severity_rank")

analysis_ready <- analysis_ready %>%
  left_join(cluster_severity %>% dplyr::select(traj_cluster_ra, severity_rank), by = "traj_cluster_ra")

# ================================================================================================
# 13) Models
# ================================================================================================

# Main site-level models using local clusters
m_traj_ra <- nnet::multinom(
  traj_cluster_ra ~ pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready
)

m_mort_ra <- glm(
  death_or_hospice ~ pm25_5y_z + no2_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_int_pm <- glm(
  death_or_hospice ~ pm25_5y_z * traj_cluster_ra + no2_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_int_no2 <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster_ra + pm25_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_los_ra <- glm(
  icu_los_hours ~ pm25_5y_z + no2_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = quasipoisson(link = "log")
)

m_severity <- MASS::polr(
  severity_rank ~ pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  Hess = TRUE
)

# landmark survival after 72h
landmark72 <- analysis_ready %>%
  mutate(
    t0_72 = admission_dttm + dhours(72),
    surv_time_ref = dplyr::coalesce(death_ts, discharge_dttm),
    eligible72 = !is.na(discharge_dttm) & discharge_dttm > t0_72,
    event_after72 = death_or_hospice & surv_time_ref > t0_72,
    time_from72_h = as.numeric(difftime(surv_time_ref, t0_72, units = "hours"))
  ) %>%
  filter(eligible72, !is.na(traj_cluster_ra), !is.na(time_from72_h), time_from72_h >= 0)

cox72 <- coxph(
  Surv(time_from72_h, event_after72) ~ traj_cluster_ra +
    pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = landmark72
)

cox72_int <- coxph(
  Surv(time_from72_h, event_after72) ~ traj_cluster_ra * no2_5y_z + pm25_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = landmark72
)

# mediation-style attenuation analyses
m_total_no2 <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_direct_no2 <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z + severity_rank +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_direct_no2_cluster <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

# standardized model output tables
model_outputs <- bind_rows(
  tidy_multinom_model(m_traj_ra, "multinom_traj_cluster"),
  tidy_glm_model(m_mort_ra, "glm_mortality"),
  tidy_glm_model(m_int_pm, "glm_interaction_pm25"),
  tidy_glm_model(m_int_no2, "glm_interaction_no2"),
  tidy_glm_model(m_los_ra, "glm_los_quasipoisson"),
  tidy_polr_model(m_severity, "polr_severity_rank"),
  tidy_glm_model(m_total_no2, "glm_total_no2"),
  tidy_glm_model(m_direct_no2, "glm_direct_no2_severity"),
  tidy_glm_model(m_direct_no2_cluster, "glm_direct_no2_cluster"),
  tidy_cox_model(cox72, "cox72_landmark"),
  tidy_cox_model(cox72_int, "cox72_landmark_interaction_no2")
)

save_model_table(model_outputs, "model_outputs_standardized")

# site-level model metadata / counts
safe_n <- function(model) {
  out <- tryCatch(nrow(model.frame(model)), error = function(e) NA_integer_)
  as.integer(out)
}

model_meta <- tibble(
  model = c(
    "multinom_traj_cluster",
    "glm_mortality",
    "glm_interaction_pm25",
    "glm_interaction_no2",
    "glm_los_quasipoisson",
    "polr_severity_rank",
    "cox72_landmark",
    "cox72_landmark_interaction_no2"
  ),
  n = c(
    safe_n(m_traj_ra),
    safe_n(m_mort_ra),
    safe_n(m_int_pm),
    safe_n(m_int_no2),
    safe_n(m_los_ra),
    safe_n(m_severity),
    nrow(landmark72),
    nrow(landmark72)
  ),
  events = c(
    NA_real_,
    sum(model.frame(m_mort_ra)$death_or_hospice, na.rm = TRUE),
    sum(model.frame(m_int_pm)$death_or_hospice, na.rm = TRUE),
    sum(model.frame(m_int_no2)$death_or_hospice, na.rm = TRUE),
    NA_real_,
    NA_real_,
    sum(landmark72$event_after72, na.rm = TRUE),
    sum(landmark72$event_after72, na.rm = TRUE)
  )
)

save_csv(model_meta, "model_metadata")

# VIF
if (requireNamespace("car", quietly = TRUE)) {
  vif_tbl <- tibble(
    term = names(car::vif(m_mort_ra)),
    vif = as.numeric(car::vif(m_mort_ra))
  )
  save_csv(vif_tbl, "model_vif_mortality")
}

# mediation-style decomposition tables
beta_total_no2  <- coef(m_total_no2)[["no2_5y_z"]]
beta_direct_no2 <- coef(m_direct_no2)[["no2_5y_z"]]

mediation_decomp <- tibble(
  exposure = "NO2_5y_z",
  beta_total = beta_total_no2,
  beta_direct = beta_direct_no2,
  beta_indirect_approx = beta_total_no2 - beta_direct_no2,
  pct_mediated_approx = (beta_total_no2 - beta_direct_no2) / beta_total_no2
)
save_csv(mediation_decomp, "mediation_style_decomposition_no2")

beta_total_no2_cluster  <- coef(m_total_no2)[["no2_5y_z"]]
beta_direct_no2_cluster <- coef(m_direct_no2_cluster)[["no2_5y_z"]]

mediation_decomp_cluster <- tibble(
  exposure = "NO2_5y_z",
  beta_total = beta_total_no2_cluster,
  beta_direct = beta_direct_no2_cluster,
  beta_indirect_approx = beta_total_no2_cluster - beta_direct_no2_cluster,
  pct_mediated_approx = (beta_total_no2_cluster - beta_direct_no2_cluster) / beta_total_no2_cluster
)
save_csv(mediation_decomp_cluster, "mediation_style_decomposition_no2_cluster")

# ================================================================================================
# 14) Predicted mortality / interaction figures
# ================================================================================================

ref_age <- mean(analysis_ready$age_years, na.rm = TRUE)
ref_charlson <- mean(analysis_ready$charlson_score, na.rm = TRUE)
ref_sofa <- mean(analysis_ready$sofa_total, na.rm = TRUE)
ref_year <- round(mean(analysis_ready$admit_year, na.rm = TRUE))
ref_sex <- analysis_ready %>% count(sex_category, sort = TRUE) %>% slice(1) %>% pull(sex_category)
ref_race <- analysis_ready %>% count(race_category, sort = TRUE) %>% slice(1) %>% pull(race_category)
ref_adv <- FALSE

# PM2.5 interaction predictions
newdat_pm <- tidyr::expand_grid(
  traj_cluster_ra = levels(droplevels(analysis_ready$traj_cluster_ra)),
  pm25_5y_z = seq(-2, 2, by = 0.1)
) %>%
  mutate(
    no2_5y_z = 0,
    age_years = ref_age,
    sex_category = ref_sex,
    race_category = ref_race,
    admit_year = ref_year,
    charlson_score = ref_charlson,
    sofa_total = ref_sofa,
    advanced_cancer_any_poa = ref_adv
  ) %>%
  mutate(pred_prob = predict(m_int_pm, newdata = ., type = "response"))

p_pm <- ggplot(newdat_pm, aes(x = pm25_5y_z, y = pred_prob, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Predicted mortality by trajectory cluster and PM2.5 exposure",
    x = "PM2.5 (5-year exposure, z-score)",
    y = "Predicted probability of in-hospital death/hospice",
    color = "Trajectory cluster"
  ) +
  theme_pub()
save_plot(p_pm, "predicted_mortality_by_cluster_pm25", w = 10, h = 7)

# NO2 interaction predictions
newdat_no2 <- tidyr::expand_grid(
  traj_cluster_ra = levels(droplevels(analysis_ready$traj_cluster_ra)),
  no2_5y_z = seq(-2, 2, by = 0.1)
) %>%
  mutate(
    pm25_5y_z = 0,
    age_years = ref_age,
    sex_category = ref_sex,
    race_category = ref_race,
    admit_year = ref_year,
    charlson_score = ref_charlson,
    sofa_total = ref_sofa,
    advanced_cancer_any_poa = ref_adv
  ) %>%
  mutate(pred_prob = predict(m_int_no2, newdata = ., type = "response"))

p_no2 <- ggplot(newdat_no2, aes(x = no2_5y_z, y = pred_prob, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Predicted mortality by trajectory cluster and NO2 exposure",
    x = "NO2 (5-year exposure, z-score)",
    y = "Predicted probability of in-hospital death/hospice",
    color = "Trajectory cluster"
  ) +
  theme_pub()
save_plot(p_no2, "predicted_mortality_by_cluster_no2", w = 10, h = 7)

contrast_pm <- newdat_pm %>%
  filter(pm25_5y_z %in% c(-1, 1)) %>%
  mutate(level = ifelse(pm25_5y_z == -1, "low", "high")) %>%
  dplyr::select(traj_cluster_ra, level, pred_prob) %>%
  pivot_wider(names_from = level, values_from = pred_prob) %>%
  mutate(
    risk_diff = high - low,
    risk_ratio = high / low
  )
save_csv(contrast_pm, "contrast_pm25_low_vs_high_by_cluster")

contrast_no2 <- newdat_no2 %>%
  filter(no2_5y_z %in% c(-1, 1)) %>%
  mutate(level = ifelse(no2_5y_z == -1, "low", "high")) %>%
  dplyr::select(traj_cluster_ra, level, pred_prob) %>%
  pivot_wider(names_from = level, values_from = pred_prob) %>%
  mutate(
    risk_diff = high - low,
    risk_ratio = high / low
  )
save_csv(contrast_no2, "contrast_no2_low_vs_high_by_cluster")

p_diff_pm <- ggplot(contrast_pm, aes(x = traj_cluster_ra, y = risk_diff, fill = traj_cluster_ra)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Absolute mortality increase from low to high PM2.5 by trajectory cluster",
    x = "Trajectory cluster",
    y = "Risk difference",
    fill = "Trajectory cluster"
  ) +
  theme_pub() +
  theme(legend.position = "none")
save_plot(p_diff_pm, "risk_difference_pm25_by_cluster", w = 9, h = 6)

p_diff_no2 <- ggplot(contrast_no2, aes(x = traj_cluster_ra, y = risk_diff, fill = traj_cluster_ra)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Absolute mortality increase from low to high NO2 by trajectory cluster",
    x = "Trajectory cluster",
    y = "Risk difference",
    fill = "Trajectory cluster"
  ) +
  theme_pub() +
  theme(legend.position = "none")
save_plot(p_diff_no2, "risk_difference_no2_by_cluster", w = 9, h = 6)

heat_pm <- newdat_pm %>% mutate(pm_bin = round(pm25_5y_z, 1))
p_heat_pm <- ggplot(heat_pm, aes(x = pm_bin, y = traj_cluster_ra, fill = pred_prob)) +
  geom_tile() +
  scale_fill_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Predicted mortality heatmap: PM2.5 × cluster",
    x = "PM2.5 (5-year exposure, z-score)",
    y = "Trajectory cluster",
    fill = "Pred. death"
  ) +
  theme_pub()
save_plot(p_heat_pm, "interaction_pm25_heatmap", w = 10, h = 5.5)

analysis_ready <- analysis_ready %>%
  mutate(no2_q = ntile(no2_5y, 4))

p_sev_no2 <- ggplot(analysis_ready, aes(x = factor(no2_q), fill = severity_rank)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Trajectory severity distribution across NO2 quartiles",
    x = "NO2 quartile",
    y = "Proportion",
    fill = "Severity rank"
  ) +
  theme_pub()
save_plot(p_sev_no2, "severity_distribution_by_no2_quartile", w = 9, h = 6)

# ================================================================================================
# 15) Additional QC / interpretation tables
# ================================================================================================

cluster_exposure_mortality <- analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    mortality = mean(death_or_hospice, na.rm = TRUE),
    pm25 = mean(pm25_5y, na.rm = TRUE),
    no2 = mean(no2_5y, na.rm = TRUE),
    .groups = "drop"
  )
save_csv(cluster_exposure_mortality, "cluster_exposure_mortality_summary")

gcs_cluster <- scores_df %>%
  left_join(analysis_ready %>% dplyr::select(hospitalization_id, traj_cluster_ra), by = "hospitalization_id") %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    median_gcs = median(numerical_value, na.rm = TRUE),
    mean_gcs = mean(numerical_value, na.rm = TRUE),
    pct_gcs_lt8 = mean(numerical_value < 8, na.rm = TRUE),
    pct_gcs_lt13 = mean(numerical_value < 13, na.rm = TRUE),
    .groups = "drop"
  )
save_csv(gcs_cluster, "cluster_gcs_summary")

neuro_prefix <- c("I61", "I62", "I63", "I64", "G93", "C793")
neuro_flags <- hospital_dx %>%
  filter(poa_present == 1) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    code = norm_code(diagnosis_code)
  ) %>%
  filter(code_matches_any_prefix(code, neuro_prefix)) %>%
  distinct(hospitalization_id) %>%
  mutate(neuro_dx = TRUE)

cluster_neuro <- analysis_ready %>%
  left_join(neuro_flags, by = "hospitalization_id") %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    neuro_dx_pct = mean(neuro_dx, na.rm = TRUE),
    .groups = "drop"
  )
save_csv(cluster_neuro, "cluster_neuro_summary")

cor_tbl <- tibble(
  variable_x = c("pm25_5y", "no2_5y"),
  variable_y = c("sofa_total", "sofa_total"),
  estimate = c(
    suppressWarnings(cor(analysis_ready$pm25_5y, analysis_ready$sofa_total, use = "complete.obs")),
    suppressWarnings(cor(analysis_ready$no2_5y, analysis_ready$sofa_total, use = "complete.obs"))
  )
)
save_csv(cor_tbl, "correlation_exposure_sofa")

# ================================================================================================
# 16) Final federated export manifest
# ================================================================================================

export_manifest <- tibble(
  category = c(
    "clustering", "clustering",
    "static_summary", "static_summary", "static_summary",
    "trajectory_summary", "trajectory_summary",
    "model_output", "model_output", "model_output",
    "figure", "figure", "figure", "figure", "figure", "figure"
  ),
  file_stub = c(
    "cluster_silhouette_ra",
    "cluster_centroids_federated",
    "cluster_static_summary",
    "cluster_sofa_summary",
    "cluster_charlson_summary",
    "cluster_signature_rs_hourly",
    "cluster_signature_arf_hourly",
    "model_outputs_standardized",
    "model_metadata",
    "mediation_style_decomposition_no2",
    "traj_ra_seqrplot_discrete",
    "cluster_sofa_domain_signature",
    "sig_rs_area_by_cluster",
    "sig_arf_area_by_cluster",
    "predicted_mortality_by_cluster_no2",
    "risk_difference_no2_by_cluster"
  ),
  description = c(
    "Silhouette diagnostics for local cluster selection.",
    "Cluster-level centroids for central harmonization.",
    "Static cluster summaries including vasoactive mean use.",
    "SOFA summaries by cluster.",
    "Charlson summaries by cluster.",
    "Hourly respiratory support proportions by cluster.",
    "Hourly ARF subtype proportions by cluster.",
    "Standardized coefficient table for all models.",
    "Model Ns and event counts.",
    "Approximate mediation-style attenuation for NO2.",
    "Representative sequence figure.",
    "SOFA domain signature figure.",
    "Respiratory support signature figure.",
    "ARF signature figure.",
    "Predicted mortality by NO2 and cluster figure.",
    "Absolute NO2 risk difference by cluster figure."
  )
)
save_csv(export_manifest, "federated_export_manifest")


# ================================================================================================
# 17) Spatial distribution of trajectory clusters (county-level federated export)
# ================================================================================================

options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

# ---------------------------
# Parameters
# ---------------------------
MIN_CELL_N <- 11L   # suppress small cells for federated export if needed
MAP_YEAR <- 2022L

# ---------------------------
# Build county-cluster table
# ---------------------------
county_cluster <- analysis_ready %>%
  mutate(
    county_code = gsub("[^0-9]", "", as.character(county_code)),
    county_code = ifelse(nchar(county_code) == 4, paste0("0", county_code), county_code)
  ) %>%
  filter(!is.na(county_code), nchar(county_code) == 5, !is.na(traj_cluster_ra)) %>%
  count(county_code, traj_cluster_ra, name = "n") %>%
  group_by(county_code) %>%
  mutate(
    county_total = sum(n),
    prop_cluster_in_county = n / county_total
  ) %>%
  ungroup() %>%
  group_by(traj_cluster_ra) %>%
  mutate(
    cluster_total = sum(n),
    prop_county_in_cluster = n / cluster_total
  ) %>%
  ungroup() %>%
  mutate(
    n_suppressed = ifelse(n < MIN_CELL_N, NA_integer_, n),
    prop_cluster_in_county_suppressed = ifelse(n < MIN_CELL_N, NA_real_, prop_cluster_in_county),
    prop_county_in_cluster_suppressed = ifelse(n < MIN_CELL_N, NA_real_, prop_county_in_cluster),
    site_name = SITE_NAME
  )

save_csv(county_cluster, "spatial_county_cluster_distribution")

# ---------------------------
# County-level summary table
# ---------------------------
county_summary <- county_cluster %>%
  group_by(county_code) %>%
  summarize(
    site_name = first(site_name),
    county_total = first(county_total),
    dominant_cluster = traj_cluster_ra[which.max(prop_cluster_in_county)],
    dominant_cluster_prop = max(prop_cluster_in_county, na.rm = TRUE),
    n_clusters_present = sum(n > 0, na.rm = TRUE),
    cluster_entropy = {
      p <- prop_cluster_in_county[prop_cluster_in_county > 0]
      -sum(p * log(p))
    },
    .groups = "drop"
  ) %>%
  mutate(
    county_total_suppressed = ifelse(county_total < MIN_CELL_N, NA_integer_, county_total),
    dominant_cluster_prop_suppressed = ifelse(county_total < MIN_CELL_N, NA_real_, dominant_cluster_prop),
    dominant_cluster_suppressed = ifelse(county_total < MIN_CELL_N, NA_character_, as.character(dominant_cluster)),
    cluster_entropy_suppressed = ifelse(county_total < MIN_CELL_N, NA_real_, cluster_entropy)
  )

save_csv(county_summary, "spatial_county_cluster_summary")

# ================================================================================================
# Maps
# ================================================================================================

# ================================================================================================
# Restrict counties to CONUS
# ================================================================================================

conus_states <- c(
  "01","04","05","06","08","09","10","11","12","13",
  "16","17","18","19","20","21","22","23","24","25",
  "26","27","28","29","30","31","32","33","34","35",
  "36","37","38","39","40","41","42","44","45","46",
  "47","48","49","50","51","53","54","55","56"
)

us_counties <- tigris::counties(cb = TRUE, year = MAP_YEAR, class = "sf") %>%
  st_transform(5070) %>%
  mutate(
    county_code = GEOID,
    state_fips = substr(GEOID, 1, 2)
  ) %>%
  filter(state_fips %in% conus_states)

# Join county summary
county_map_df <- us_counties %>%
  left_join(county_summary, by = "county_code")

# ---------------------------
# Map 1: dominant cluster by county
# ---------------------------
p_map_dominant <- ggplot(county_map_df) +
  geom_sf(aes(fill = dominant_cluster_suppressed), color = NA) +
  labs(
    title = "Dominant respiratory trajectory cluster by county",
    subtitle = paste0("Site: ", SITE_NAME),
    fill = "Dominant cluster"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot(p_map_dominant, "map_county_dominant_cluster", w = 12, h = 8)

# ---------------------------
# Map 2: dominant cluster concentration
# ---------------------------
p_map_concentration <- ggplot(county_map_df) +
  geom_sf(aes(fill = dominant_cluster_prop_suppressed), color = NA) +
  scale_fill_viridis_c(na.value = "grey90") +
  labs(
    title = "Dominant cluster concentration by county",
    subtitle = paste0("Site: ", SITE_NAME),
    fill = "Dominant cluster proportion"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot(p_map_concentration, "map_county_dominant_cluster_prop", w = 12, h = 8)

# ---------------------------
# Map 3: cluster heterogeneity (entropy)
# ---------------------------
p_map_entropy <- ggplot(county_map_df) +
  geom_sf(aes(fill = cluster_entropy_suppressed), color = NA) +
  scale_fill_viridis_c(na.value = "grey90") +
  labs(
    title = "Trajectory cluster heterogeneity by county",
    subtitle = paste0("Site: ", SITE_NAME),
    fill = "Entropy"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot(p_map_entropy, "map_county_cluster_entropy", w = 12, h = 8)

# ================================================================================================
# Optional: one map per cluster showing within-county prevalence
# ================================================================================================

cluster_map_long <- county_cluster %>%
  transmute(
    site_name,
    county_code,
    traj_cluster_ra,
    county_total,
    n = n_suppressed,
    prop_cluster_in_county = prop_cluster_in_county_suppressed
  ) %>%
  left_join(us_counties %>% dplyr::select(county_code, geometry), by = "county_code") %>%
  st_as_sf()

p_map_cluster_facets <- ggplot(cluster_map_long) +
  geom_sf(aes(fill = prop_cluster_in_county), color = NA) +
  scale_fill_viridis_c(na.value = "grey90") +
  facet_wrap(~ traj_cluster_ra) +
  labs(
    title = "Within-county prevalence of trajectory clusters",
    subtitle = paste0("Site: ", SITE_NAME),
    fill = "Cluster proportion"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot(p_map_cluster_facets, "map_county_cluster_prevalence_facets", w = 14, h = 10)

# ================================================================================================
# 18) CONSORT-style cohort flow table
# ================================================================================================

# Uses objects created in 01_lungcx_cohort.R:
#   flow_lung
#   exclusion_lung
#   cohort_lung

consort_flow_counts <- flow_lung %>%
  transmute(
    site_name = SITE_NAME,
    step,
    step_order = row_number(),
    n_remaining = remaining,
    n_excluded_at_step = excluded_at_step
  )

save_csv(consort_flow_counts, "consort_flow_counts")

consort_exclusion_reasons <- exclusion_lung %>%
  count(reason, name = "n") %>%
  mutate(site_name = SITE_NAME) %>%
  relocate(site_name)

save_csv(consort_exclusion_reasons, "consort_exclusion_reasons")

# ================================================================================================
# 19) Table 1 exports (federated-friendly long format)
# ================================================================================================

# analysis_ready already contains from 01_lungcx_cohort.R:
#   hospitalization_id
#   admission_dttm, discharge_dttm
#   death_in_hosp, hospice_discharge, death_or_hospice
#   icu_los_hours
#   age_years, sex_category, race_category, ethnicity_category
#   census_tract, county_code
#   device_t0
#   imv_hours_icu, imv_hours_72h
#   any_imv_icu, any_imv_72h
#   escalations_icu, deescalations_icu
#   pm25_1y, pm25_3y, pm25_5y
#   no2_1y, no2_3y, no2_5y
#
# and from 02_federated_clusters.R gets augmented with:
#   traj_cluster_ra
#   vaso_any_72h, vaso_hours_72h, vaso_mean_72h
#   charlson_score
#   advanced_cancer_any_poa, metastatic_cancer_poa,
#   malignant_pleural_effusion_poa, cachexia_poa
#   sofa_total, sofa_resp, sofa_cv, sofa_renal, sofa_liver, sofa_coag, sofa_cns
#   severity_rank

table1_dat <- analysis_ready %>%
  left_join(vaso_static, by = "hospitalization_id") %>%
  mutate(
    site_name = SITE_NAME,
    
    age_years = as.numeric(age_years),
    icu_los_hours = as.numeric(icu_los_hours),
    
    imv_hours_icu = as.numeric(imv_hours_icu),
    imv_hours_72h = as.numeric(imv_hours_72h),
    escalations_icu = as.numeric(escalations_icu),
    deescalations_icu = as.numeric(deescalations_icu),
    
    pm25_1y = as.numeric(pm25_1y),
    pm25_3y = as.numeric(pm25_3y),
    pm25_5y = as.numeric(pm25_5y),
    no2_1y  = as.numeric(no2_1y),
    no2_3y  = as.numeric(no2_3y),
    no2_5y  = as.numeric(no2_5y),
    
    charlson_score = as.numeric(charlson_score),
    sofa_total = as.numeric(sofa_total),
    sofa_resp = as.numeric(sofa_resp),
    sofa_cv = as.numeric(sofa_cv),
    sofa_renal = as.numeric(sofa_renal),
    sofa_liver = as.numeric(sofa_liver),
    sofa_coag = as.numeric(sofa_coag),
    sofa_cns = as.numeric(sofa_cns),
    
    death_in_hosp = as.integer(coalesce(death_in_hosp, FALSE)),
    hospice_discharge = as.integer(coalesce(hospice_discharge, FALSE)),
    death_or_hospice = as.integer(coalesce(death_or_hospice, FALSE)),
    any_imv_icu = as.integer(coalesce(any_imv_icu, FALSE)),
    any_imv_72h = as.integer(coalesce(any_imv_72h, FALSE)),
    vaso_any_72h = as.integer(coalesce(vaso_any_72h, 0L)),
    
    advanced_cancer_any_poa = as.integer(coalesce(advanced_cancer_any_poa, FALSE)),
    metastatic_cancer_poa = as.integer(coalesce(metastatic_cancer_poa, FALSE)),
    malignant_pleural_effusion_poa = as.integer(coalesce(malignant_pleural_effusion_poa, FALSE)),
    cachexia_poa = as.integer(coalesce(cachexia_poa, FALSE))
  )

# ---------------------------
# Summary helpers
# ---------------------------
summ_cont <- function(data, var, by) {
  v <- rlang::sym(var)
  b <- rlang::sym(by)
  
  data %>%
    group_by(!!b) %>%
    summarize(
      n_nonmissing = sum(!is.na(!!v)),
      mean = mean(!!v, na.rm = TRUE),
      sd = sd(!!v, na.rm = TRUE),
      median = median(!!v, na.rm = TRUE),
      p25 = quantile(!!v, 0.25, na.rm = TRUE),
      p75 = quantile(!!v, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      site_name = SITE_NAME,
      variable = var,
      level = NA_character_,
      summary_type = "continuous"
    ) %>%
    rename(stratum = !!b) %>%
    relocate(site_name, stratum, variable, level, summary_type)
}

summ_cat <- function(data, var, by) {
  v <- rlang::sym(var)
  b <- rlang::sym(by)
  
  data %>%
    mutate(.level = as.character(!!v)) %>%
    group_by(!!b, .level) %>%
    summarize(n = n(), .groups = "drop_last") %>%
    mutate(
      denom = sum(n),
      prop = n / denom
    ) %>%
    ungroup() %>%
    mutate(
      site_name = SITE_NAME,
      variable = var,
      level = .level,
      summary_type = "categorical"
    ) %>%
    rename(stratum = !!b) %>%
    dplyr::select(site_name, stratum, variable, level, summary_type, n, denom, prop)
}

summ_cont_overall <- function(data, var) {
  v <- rlang::sym(var)
  
  data %>%
    summarize(
      n_nonmissing = sum(!is.na(!!v)),
      mean = mean(!!v, na.rm = TRUE),
      sd = sd(!!v, na.rm = TRUE),
      median = median(!!v, na.rm = TRUE),
      p25 = quantile(!!v, 0.25, na.rm = TRUE),
      p75 = quantile(!!v, 0.75, na.rm = TRUE)
    ) %>%
    mutate(
      site_name = SITE_NAME,
      stratum = "Overall",
      variable = var,
      level = NA_character_,
      summary_type = "continuous"
    ) %>%
    relocate(site_name, stratum, variable, level, summary_type)
}

summ_cat_overall <- function(data, var) {
  v <- rlang::sym(var)
  
  data %>%
    mutate(.level = as.character(!!v)) %>%
    count(.level, name = "n") %>%
    mutate(
      denom = sum(n),
      prop = n / denom,
      site_name = SITE_NAME,
      stratum = "Overall",
      variable = var,
      level = .level,
      summary_type = "categorical"
    ) %>%
    dplyr::select(site_name, stratum, variable, level, summary_type, n, denom, prop)
}

# ---------------------------
# Variables for Table 1
# ---------------------------
cont_vars <- c(
  "age_years",
  "icu_los_hours",
  "imv_hours_icu",
  "imv_hours_72h",
  "escalations_icu",
  "deescalations_icu",
  "vaso_hours_72h",
  "vaso_mean_72h",
  "pm25_1y",
  "pm25_3y",
  "pm25_5y",
  "no2_1y",
  "no2_3y",
  "no2_5y",
  "charlson_score",
  "sofa_total",
  "sofa_resp",
  "sofa_cv",
  "sofa_renal",
  "sofa_liver",
  "sofa_coag",
  "sofa_cns"
)

cat_vars <- c(
  "sex_category",
  "race_category",
  "ethnicity_category",
  "device_t0",
  "traj_cluster_ra",
  "severity_rank",
  "death_in_hosp",
  "hospice_discharge",
  "death_or_hospice",
  "any_imv_icu",
  "any_imv_72h",
  "vaso_any_72h",
  "advanced_cancer_any_poa",
  "metastatic_cancer_poa",
  "malignant_pleural_effusion_poa",
  "cachexia_poa"
)

# ---------------------------
# Overall Table 1
# ---------------------------
table1_overall_cont <- purrr::map_dfr(cont_vars, ~ summ_cont_overall(table1_dat, .x))
table1_overall_cat  <- purrr::map_dfr(cat_vars, ~ summ_cat_overall(table1_dat, .x))

save_csv(table1_overall_cont, "table1_overall_continuous")
save_csv(table1_overall_cat, "table1_overall_categorical")

# ---------------------------
# Table 1 by trajectory cluster
# ---------------------------
table1_by_cluster_cont <- purrr::map_dfr(cont_vars, ~ summ_cont(table1_dat, .x, by = "traj_cluster_ra"))
table1_by_cluster_cat  <- purrr::map_dfr(cat_vars, ~ summ_cat(table1_dat, .x, by = "traj_cluster_ra"))

save_csv(table1_by_cluster_cont, "table1_by_cluster_continuous")
save_csv(table1_by_cluster_cat, "table1_by_cluster_categorical")

# ---------------------------
# Table 1 by NO2 quartile
# ---------------------------
table1_dat <- table1_dat %>%
  mutate(
    no2_quartile = ntile(no2_5y, 4),
    no2_quartile = factor(no2_quartile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4"))
  )

table1_by_no2_cont <- purrr::map_dfr(cont_vars, ~ summ_cont(table1_dat, .x, by = "no2_quartile"))
table1_by_no2_cat  <- purrr::map_dfr(cat_vars, ~ summ_cat(table1_dat, .x, by = "no2_quartile"))

save_csv(table1_by_no2_cont, "table1_by_no2_continuous")
save_csv(table1_by_no2_cat, "table1_by_no2_categorical")



cat("\nDone. All trajectory panels, clustering outputs, standardized model tables, and federated summaries were saved to:\n",
    out_dir, "\n", sep = "")
cat("\nPlease upload the output folder to Box.\n")
