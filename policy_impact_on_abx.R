##### Policy impact on unnecessary antibiotics #####
# Assumes create_new_pros_model.R and sbi_decision_policy_eplx.R have already been run
# and objects below already exist in the current R session:
#   - pros_one_model
#   - test_decisions_best

library(data.table)
library(dplyr)
library(lubridate)

### Create flag for whether patient received any antibiotics in the PICU during the first 24 hours.
abx_raw <- read.csv(file = "/phi/sbi/sbi_blake/antinfective_export_pros_091225.csv")

abx_df <- abx_raw %>%
  mutate(
    order_inst = as.POSIXct(order_inst, tz = "UTC"),
    taken_time = as.POSIXct(taken_time, tz = "UTC")
  )

abx_df <- abx_df %>%
  mutate(
    order_inst = force_tz(order_inst, tzone = "America/Denver"),
    taken_time = force_tz(taken_time, tzone = "America/Denver")
  )

abx_df$PAT_MRN_ID <- as.character(abx_df$PAT_MRN_ID)

abx_df <- abx_df %>%
  rename(
    mrn = PAT_MRN_ID,
    hsp_account_id = HSP_ACCOUNT_ID
  )

colnames(abx_df) <- str_to_lower(colnames(abx_df))

setDT(abx_df)
abx_df[, taken_time := as.POSIXct(taken_time)]

# ------------------------------------------------------------
# Create flag for any antibiotic dose during first 24h of PICU
# ------------------------------------------------------------

# Make a one-row-per-PICU-encounter dataset
picu_24h_windows <- pros_one_model %>%
  dplyr::distinct(
    study_id,
    pat_enc_csn_id,
    picu_adm_date_time
  ) %>%
  dplyr::filter(
    !is.na(study_id),
    !is.na(pat_enc_csn_id),
    !is.na(picu_adm_date_time)
  )

# Convert to data.table
data.table::setDT(picu_24h_windows)

# Ensure PICU admission time is POSIXct
picu_24h_windows[, picu_adm_date_time := as.POSIXct(picu_adm_date_time)]

# Create the first-24h PICU window
picu_24h_windows[, window_start := picu_adm_date_time]
picu_24h_windows[, window_end   := picu_adm_date_time + lubridate::hours(24)]

# Create antibiotic point-event table
abx_events <- data.table::copy(abx_df)

data.table::setDT(abx_events)

abx_events[, taken_time := as.POSIXct(taken_time)]

abx_events <- abx_events[
  !is.na(pat_enc_csn_id) &
    !is.na(taken_time)
]

# Treat each antibiotic administration as a point event
abx_events[, abx_start := taken_time]
abx_events[, abx_end   := taken_time]

# Key the PICU-window table
data.table::setkey(
  picu_24h_windows,
  pat_enc_csn_id,
  window_start,
  window_end
)

# Key the antibiotic-event table
data.table::setkey(
  abx_events,
  pat_enc_csn_id,
  abx_start,
  abx_end
)

# Find antibiotic administrations that occurred within first 24h of PICU admission
abx_first_24h_overlaps <- data.table::foverlaps(
  x = abx_events,
  y = picu_24h_windows,
  by.x = c("pat_enc_csn_id", "abx_start", "abx_end"),
  by.y = c("pat_enc_csn_id", "window_start", "window_end"),
  type = "within",
  nomatch = 0L
)

# Identify first antibiotic dose time within first 24h PICU window per encounter
first_abx_in_first_24h <- abx_first_24h_overlaps[, .(
  first_abx_time_in_first_24h_picu = min(taken_time, na.rm = TRUE)
), by = .(study_id, pat_enc_csn_id, picu_adm_date_time)]

first_abx_in_first_24h[, first_abx_hour_in_first_24h_picu :=
  as.numeric(difftime(first_abx_time_in_first_24h_picu, picu_adm_date_time, units = "hours"))
]

# Create one-row-per-study_id flag
abx_first_24h_flag <- picu_24h_windows[, .(
  study_id,
  pat_enc_csn_id,
  picu_adm_date_time
)]

abx_first_24h_flag[, abx_in_first_24h_picu := 0L]

study_ids_with_abx_first_24h <- unique(abx_first_24h_overlaps$study_id)

abx_first_24h_flag[
  study_id %in% study_ids_with_abx_first_24h,
  abx_in_first_24h_picu := 1L
]

# Join flag back onto pros_one_model
pros_one_model <- pros_one_model %>%
  dplyr::left_join(
    abx_first_24h_flag %>%
      dplyr::as_tibble() %>%
      dplyr::select(
        study_id,
        abx_in_first_24h_picu
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    first_abx_in_first_24h %>%
      dplyr::as_tibble() %>%
      dplyr::select(
        study_id,
        first_abx_time_in_first_24h_picu,
        first_abx_hour_in_first_24h_picu
      ),
    by = "study_id"
  ) %>%
  dplyr::mutate(
    abx_in_first_24h_picu = dplyr::if_else(
      is.na(abx_in_first_24h_picu),
      0L,
      abx_in_first_24h_picu
    )
  )




library(dplyr)
library(tibble)

cap_days <- 7

# Keep only policy-predicted SBI-negative encounters on the prospective test set
policy_sbi_negative <- test_decisions_best %>%
  dplyr::filter(final_state == "ruled_out") %>%
  dplyr::mutate(study_id = as.character(study_id))

# Pull antibiotic duration at the exact time the policy first rules out SBI
# (decision_hour), then cap durations at 7 days.
abx_at_decision <- pros_one_model %>%
  dplyr::transmute(
    study_id = as.character(study_id),
    hours_since_picu_adm = as.integer(hours_since_picu_adm),
    abx_duration_after_score = as.numeric(abx_duration_after_score)
  ) %>%
  dplyr::group_by(study_id, hours_since_picu_adm) %>%
  dplyr::summarise(
    abx_duration_after_score = dplyr::first(abx_duration_after_score),
    .groups = "drop"
  )

policy_abx <- policy_sbi_negative %>%
  dplyr::left_join(
    abx_at_decision,
    by = c("study_id", "decision_hour" = "hours_since_picu_adm")
  ) %>%
  dplyr::mutate(
    abx_duration_after_score_capped = dplyr::if_else(
      is.na(abx_duration_after_score),
      NA_real_,
      pmin(abx_duration_after_score, cap_days)
    )
  )

# Restrict "preventable unnecessary antibiotics" numerator calculations to true
# SBI-negative encounters that the policy ruled out.
policy_true_neg <- policy_abx %>%
  dplyr::filter(sbi_present == 0)

# Denominator cohort: all true SBI-negative encounters (regardless of final policy
# state), using each encounter's score timepoint (decision_hour), with any
# antibiotic exposure after that score timepoint.
sbi_negative_all_abx <- test_decisions_best %>%
  dplyr::filter(sbi_present == 0) %>%
  dplyr::mutate(study_id = as.character(study_id)) %>%
  dplyr::left_join(
    abx_at_decision,
    by = c("study_id", "decision_hour" = "hours_since_picu_adm")
  ) %>%
  dplyr::mutate(
    abx_duration_after_score_capped = dplyr::if_else(
      is.na(abx_duration_after_score),
      NA_real_,
      pmin(abx_duration_after_score, cap_days)
    )
  )

# 1) Timing of SBI-negative prediction (rule-out hour)
metric_1_timing <- policy_abx %>%
  dplyr::summarise(
    n_predicted_sbi_negative = dplyr::n(),
    median_ruleout_hour = stats::median(decision_hour, na.rm = TRUE),
    q1_ruleout_hour = stats::quantile(decision_hour, probs = 0.25, na.rm = TRUE),
    q3_ruleout_hour = stats::quantile(decision_hour, probs = 0.75, na.rm = TRUE)
  )

# 2) Number and proportion with potentially preventable unnecessary antibiotics
# Numerator: true SBI-negative ruled-out encounters with >0 capped abx duration after score
# Denominator: all true SBI-negative encounters with any (>0) abx duration after score.
metric_2_preventable_count_prop <- dplyr::summarise(
  policy_true_neg,
  n_sbi_negative_with_preventable_unnecessary_abx = sum(abx_duration_after_score_capped > 0, na.rm = TRUE)
) %>%
  dplyr::bind_cols(
    sbi_negative_all_abx %>%
  dplyr::summarise(
    n_sbi_negative_with_any_abx_after_score = sum(abx_duration_after_score_capped > 0, na.rm = TRUE)
  ) %>%
      dplyr::select(n_sbi_negative_with_any_abx_after_score)
  ) %>%
  dplyr::mutate(
    prop_preventable_unnecessary_abx = dplyr::if_else(
      n_sbi_negative_with_any_abx_after_score > 0,
      n_sbi_negative_with_preventable_unnecessary_abx / n_sbi_negative_with_any_abx_after_score,
      NA_real_
    )
  )

# 3) Distribution of potentially preventable unnecessary antibiotic duration (days)
metric_3_preventable_duration <- policy_true_neg %>%
  dplyr::filter(abx_duration_after_score_capped > 0) %>%
  dplyr::summarise(
    n_with_preventable_unnecessary_abx = dplyr::n(),
    median_preventable_abx_days = stats::median(abx_duration_after_score_capped, na.rm = TRUE),
    q1_preventable_abx_days = stats::quantile(abx_duration_after_score_capped, probs = 0.25, na.rm = TRUE),
    q3_preventable_abx_days = stats::quantile(abx_duration_after_score_capped, probs = 0.75, na.rm = TRUE)
  )

# 4) Among SBI-negative encounters with antibiotics in first 24h of PICU admission,
# proportion identified by policy as ruled out before first PICU antibiotic dose.
sbi_negative_first_24h_abx <- test_decisions_best %>%
  dplyr::filter(sbi_present == 0) %>%
  dplyr::mutate(study_id = as.character(study_id)) %>%
  dplyr::left_join(
    pros_one_model %>%
      dplyr::transmute(
        study_id = as.character(study_id),
        abx_in_first_24h_picu,
        first_abx_hour_in_first_24h_picu
      ) %>%
      dplyr::distinct(study_id, .keep_all = TRUE),
    by = "study_id"
  ) %>%
  dplyr::filter(abx_in_first_24h_picu == 1, !is.na(first_abx_hour_in_first_24h_picu))

metric_4_identified_before_first_24h_abx <- sbi_negative_first_24h_abx %>%
  dplyr::summarise(
    n_sbi_negative_with_abx_in_first_24h_picu = dplyr::n(),
    n_identified_prior_to_first_picu_abx_dose = sum(final_state == "ruled_out" & decision_hour < first_abx_hour_in_first_24h_picu, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    prop_identified_prior_to_first_picu_abx_dose = dplyr::if_else(
      n_sbi_negative_with_abx_in_first_24h_picu > 0,
      n_identified_prior_to_first_picu_abx_dose / n_sbi_negative_with_abx_in_first_24h_picu,
      NA_real_
    )
  )

# Combined output table and component tables for easy downstream use
policy_impact_on_abx <- list(
  model_name = "prospective_random_forest",
  cap_days = cap_days,
  metric_1_timing = metric_1_timing,
  metric_2_preventable_count_prop = metric_2_preventable_count_prop,
  metric_3_preventable_duration = metric_3_preventable_duration,
  metric_4_identified_before_first_24h_abx = metric_4_identified_before_first_24h_abx,
  patient_level = policy_abx
)

metric_1_timing
metric_2_preventable_count_prop
metric_3_preventable_duration
metric_4_identified_before_first_24h_abx
