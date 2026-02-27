





suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(TraMineR)
  library(cluster)
})

# --------------------------- Trajectory: hourly device state (0–72h) ---------------------------


H_MAX <- 72L

keys <- cohort_lung %>%
  transmute(hospitalization_id = as.character(hospitalization_id),
            t0 = as.POSIXct(t0, tz="UTC"))

rs_for_traj <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = as.POSIXct(recorded_dttm, tz="UTC"),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h    = floor(dt_h),
    dev  = device_state_clif(device_category),
    rank = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

# Dominant device per hour = max rank; if ties, take first after arranging by recorded_dttm
hourly_state <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
  select(hospitalization_id, h, state = dev)

# Complete 0..H_MAX and fill missing with LOCF/NOCB
panel <- hourly_state %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(state, .direction = "down") %>%   # LOCF
  tidyr::fill(state, .direction = "up") %>%     # NOCB
  mutate(state = coalesce(state, "UNK")) %>%
  ungroup()

# Wide
seq_wide <- panel %>%
  mutate(h_col = paste0("H", h)) %>%
  select(hospitalization_id, h_col, state) %>%
  pivot_wider(names_from = h_col, values_from = state)

# Alphabet in clinical order
alphabet <- c("ROOM AIR","NASAL CANNULA","FACE MASK","TRACH COLLAR","HIGH FLOW NC","CPAP","NIPPV","IMV","OTHER","UNK")

seq_cols <- seq_wide %>% select(matches("^H\\d+$"))

seq_obj <- TraMineR::seqdef(
  seq_cols,
  alphabet = alphabet,
  states   = alphabet,
  labels   = alphabet,
  xtstep   = 6
)

# Distance + clustering
dist_om <- TraMineR::seqdist(seq_obj, method = "OM", indel = 1, sm = "TRATE")
hc <- hclust(as.dist(dist_om), method = "ward.D2")

K <- 5
cluster <- cutree(hc, k = K)

traj_assign <- tibble(
  hospitalization_id = seq_wide$hospitalization_id,
  traj_cluster = as.integer(cluster)
)


# Plots (base graphics)
png(file.path(out_dir, make_name("traj_seqdplot_by_cluster", "png")), width=1400, height=1400, res=160)
TraMineR::seqdplot(seq_obj, group = cluster, border = NA, main = "CLIF respiratory support distribution by cluster (0–72h)")
dev.off()

png(file.path(out_dir, make_name("traj_seqrplot_by_cluster", "png")), width=1500, height=1700, res=160)
TraMineR::seqrplot(seq_obj, group = cluster, method = "OM", indel = 1, sm = "TRATE",
                   main = "Representative CLIF respiratory support sequences by cluster")
dev.off()

# Summaries
traj_summary <- analysis_ready %>%
  left_join(traj_assign, by = "hospitalization_id") %>%
  group_by(traj_cluster) %>%
  summarize(
    n = n(),
    death_or_hospice = mean(death_or_hospice, na.rm=TRUE),
    icu_los_h_median = median(icu_los_hours, na.rm=TRUE),
    any_imv_72h = mean(any_imv_72h, na.rm=TRUE),
    imv_hours_72h_median = median(imv_hours_72h, na.rm=TRUE),
    pm25_3y_median = median(pm25_3y, na.rm=TRUE),
    no2_3y_median  = median(no2_3y,  na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(n))

save_csv(traj_summary, "traj_cluster_summary")

vitals <- get_min("vitals", c(
  "hospitalization_id","recorded_dttm","vital_category","vital_value"
)) %>%
  rename_with(tolower) %>%
  mutate(
    recorded_dttm   = safe_posix(recorded_dttm),
    vital_category  = tolower(as.character(vital_category)),
    vital_value     = as_num(vital_value)
  )

H_MAX <- 72L

keys <- cohort_lung %>%
  transmute(hospitalization_id = as.character(hospitalization_id),
            t0 = as.POSIXct(t0, tz="UTC"))

vitals_hourly <- vitals %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm,
    vital_category,
    vital_value
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h = floor(dt_h)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX) %>%
  # keep only categories we care about for the trajectory
  filter(vital_category %in% c("spo2","respiratory_rate","map","heart_rate")) %>%
  group_by(hospitalization_id, h, vital_category) %>%
  summarize(v = suppressWarnings(median(vital_value, na.rm=TRUE)), .groups="drop") %>%
  mutate(v = ifelse(is.infinite(v), NA, v)) %>%
  tidyr::pivot_wider(names_from = vital_category, values_from = v) %>%
  rename(
    spo2_h = spo2,
    rr_h   = respiratory_rate,
    map_h  = map,
    hr_h   = heart_rate
  )

# assumes device_state_clif() and device_rank_clif() are defined
rs_for_traj <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h    = floor(dt_h),
    dev  = device_state_clif(device_category),
    rank = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

hourly_state <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(hospitalization_id, h, rs_state = dev)

rs_panel <- hourly_state %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(rs_state, .direction = "down") %>%
  tidyr::fill(rs_state, .direction = "up") %>%
  mutate(rs_state = coalesce(rs_state, "UNK")) %>%
  ungroup()

phys_panel <- rs_panel %>%
  left_join(vitals_hourly, by=c("hospitalization_id","h")) %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  # keep RS fully populated
  tidyr::fill(rs_state, .direction = "downup") %>%
  ungroup()

spo2_bin <- function(x) {
  case_when(
    is.na(x)      ~ "SPO2_UNK",
    x < 88        ~ "SPO2_VLOW",
    x >= 88 & x < 92 ~ "SPO2_LOW",
    x >= 92 & x < 96 ~ "SPO2_OK",
    x >= 96       ~ "SPO2_HI"
  )
}

rr_bin <- function(x) {
  case_when(
    is.na(x)         ~ "RR_UNK",
    x < 12           ~ "RR_LOW",
    x >= 12 & x < 20 ~ "RR_NL",
    x >= 20 & x < 30 ~ "RR_HI",
    x >= 30          ~ "RR_VHI"
  )
}

map_bin <- function(x) {
  case_when(
    is.na(x)          ~ "MAP_UNK",
    x < 65            ~ "MAP_LOW",
    x >= 65 & x < 80  ~ "MAP_OK",
    x >= 80           ~ "MAP_HI"
  )
}

phys_panel <- phys_panel %>%
  mutate(
    spo2_cat = spo2_bin(spo2_h),
    rr_cat   = rr_bin(rr_h),
    map_cat  = map_bin(map_h)   # optional: keep if you want hemodynamics in the state
  )

phys_panel <- phys_panel %>%
  mutate(
    comp_state = paste(rs_state, spo2_cat, rr_cat, sep="|")
    # if adding MAP: paste(rs_state, spo2_cat, rr_cat, map_cat, sep="|")
  )

seq_wide_phys <- phys_panel %>%
  mutate(h_col = paste0("H", h)) %>%
  select(hospitalization_id, h_col, comp_state) %>%
  pivot_wider(names_from = h_col, values_from = comp_state)

seq_cols_phys <- seq_wide_phys %>% select(matches("^H\\d+$"))

alphabet_phys <- sort(unique(unlist(seq_cols_phys)))

seq_obj_phys <- TraMineR::seqdef(
  seq_cols_phys,
  alphabet = alphabet_phys,
  states   = alphabet_phys,
  labels   = alphabet_phys,
  xtstep   = 6
)

dist_om_phys <- TraMineR::seqdist(seq_obj_phys, method="OM", indel=1, sm="TRATE")
hc_phys <- hclust(as.dist(dist_om_phys), method="ward.D2")

K <- 6
cluster_phys <- cutree(hc_phys, k = K)

traj_assign_phys <- tibble(
  hospitalization_id = seq_wide_phys$hospitalization_id,
  traj_cluster_phys  = as.integer(cluster_phys)
)

save_csv(traj_assign_phys, "traj_cluster_assignments_phys")

marg_rs <- phys_panel %>%
  left_join(traj_assign_phys, by="hospitalization_id") %>%
  group_by(traj_cluster_phys, h, rs_state) %>%
  summarize(p = n() / sum(n()), .groups="drop")

marg_spo2 <- phys_panel %>%
  left_join(traj_assign_phys, by="hospitalization_id") %>%
  group_by(traj_cluster_phys, h, spo2_cat) %>%
  summarize(p = n() / sum(n()), .groups="drop")

marg_rr <- phys_panel %>%
  left_join(traj_assign_phys, by="hospitalization_id") %>%
  group_by(traj_cluster_phys, h, rr_cat) %>%
  summarize(p = n() / sum(n()), .groups="drop")

save_csv(marg_rs,   "traj_phys_marginal_rs_by_cluster_hour")
save_csv(marg_spo2, "traj_phys_marginal_spo2_by_cluster_hour")
save_csv(marg_rr,   "traj_phys_marginal_rr_by_cluster_hour")

p1 <- ggplot(marg_rs, aes(x=h, y=p, fill=rs_state)) +
  geom_area() +
  facet_wrap(~traj_cluster_phys, ncol=2) +
  labs(title="Respiratory support state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
  theme_pub()

p1

p2 <- ggplot(marg_spo2, aes(x=h, y=p, fill=spo2_cat)) +
  geom_area() +
  facet_wrap(~traj_cluster_phys, ncol=2) +
  labs(title="SpO2 state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
  theme_pub()

p2

p3 <- ggplot(marg_rr, aes(x=h, y=p, fill=rr_cat)) +
  geom_area() +
  facet_wrap(~traj_cluster_phys, ncol=2) +
  labs(title="Resp rate state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
  theme_pub()

p3

save_plot(p1, "traj_phys_rs_area", w=11, h=7)

















