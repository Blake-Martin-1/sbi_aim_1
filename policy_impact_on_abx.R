##### Policy impact on unnecessary antibiotics #####
# Assumes create_new_pros_model.R and sbi_decision_policy_eplx.R have already been run
# and objects below already exist in the current R session:
#   - pros_one_model
#   - test_decisions_best

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

# Combined output table and component tables for easy downstream use
policy_impact_on_abx <- list(
  model_name = "prospective_random_forest",
  cap_days = cap_days,
  metric_1_timing = metric_1_timing,
  metric_2_preventable_count_prop = metric_2_preventable_count_prop,
  metric_3_preventable_duration = metric_3_preventable_duration,
  patient_level = policy_abx
)

metric_1_timing
metric_2_preventable_count_prop
metric_3_preventable_duration
