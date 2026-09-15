### Create the chart-review roster for false-negative decision-policy results.

create_false_negative_chart_review_table <- function(
    decisions,
    prospective_cohort,
    expected_n = NULL) {
  required_decision_cols <- c(
    "study_id", "sbi_present", "final_state", "decision_hour"
  )
  required_cohort_cols <- c("study_id", "picu_adm_date_time")

  missing_decision_cols <- setdiff(required_decision_cols, names(decisions))
  missing_cohort_cols <- setdiff(required_cohort_cols, names(prospective_cohort))
  if (length(missing_decision_cols) > 0L) {
    stop(
      "Decision data are missing required columns: ",
      paste(missing_decision_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(missing_cohort_cols) > 0L) {
    stop(
      "Prospective cohort data are missing required columns: ",
      paste(missing_cohort_cols, collapse = ", "),
      call. = FALSE
    )
  }

  admission_times <- prospective_cohort %>%
    dplyr::transmute(
      study_id = as.character(study_id),
      picu_adm_date_time = picu_adm_date_time
    ) %>%
    dplyr::filter(!is.na(study_id), !is.na(picu_adm_date_time)) %>%
    dplyr::distinct()

  duplicate_times <- admission_times %>%
    dplyr::count(study_id, name = "n_admission_times") %>%
    dplyr::filter(n_admission_times > 1L)
  if (nrow(duplicate_times) > 0L) {
    stop(
      "Multiple PICU admission times were found for study_id(s): ",
      paste(duplicate_times$study_id, collapse = ", "),
      call. = FALSE
    )
  }

  chart_review_table <- decisions %>%
    dplyr::transmute(
      study_id = as.character(study_id),
      sbi_present = as.integer(sbi_present),
      final_state = as.character(final_state),
      hours_after_picu_admission = as.numeric(decision_hour)
    ) %>%
    dplyr::filter(sbi_present == 1L, final_state == "ruled_out") %>%
    dplyr::select(study_id, hours_after_picu_admission) %>%
    dplyr::left_join(admission_times, by = "study_id") %>%
    dplyr::select(
      study_id,
      picu_adm_date_time,
      hours_after_picu_admission
    ) %>%
    dplyr::arrange(hours_after_picu_admission, study_id)

  if (anyNA(chart_review_table$picu_adm_date_time)) {
    missing_ids <- chart_review_table$study_id[
      is.na(chart_review_table$picu_adm_date_time)
    ]
    stop(
      "PICU admission time was not found for study_id(s): ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(expected_n) && nrow(chart_review_table) != expected_n) {
    stop(
      "Expected ", expected_n, " false-negative test admissions, but found ",
      nrow(chart_review_table), ".",
      call. = FALSE
    )
  }

  chart_review_table
}

# These objects are created by create_new_pros_model.R and
# sbi_decision_policy_explore.R, which aim_1_main.R sources first.
false_negative_chart_review_table <- create_false_negative_chart_review_table(
  decisions = test_decisions_best,
  prospective_cohort = pros_one_model,
  expected_n = 18L
)

false_negative_chart_review_output_path <- getOption(
  "false_negative_chart_review_output_path",
  file.path(sbi_blake_phi_path, "false_negative_chart_review_table.csv")
)

readr::write_csv(
  false_negative_chart_review_table,
  false_negative_chart_review_output_path
)

message(
  "Wrote ", nrow(false_negative_chart_review_table),
  " false-negative admissions to ",
  false_negative_chart_review_output_path
)
