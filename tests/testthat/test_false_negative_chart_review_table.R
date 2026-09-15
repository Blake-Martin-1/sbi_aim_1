library(testthat)
library(dplyr)

script_lines <- readLines("../../create_false_negative_chart_review_table.R")
function_end <- grep("^# These objects are created", script_lines) - 1L
eval(parse(text = script_lines[seq_len(function_end)]))

test_that("chart-review table contains only false-negative admissions", {
  decisions <- tibble::tribble(
    ~study_id, ~sbi_present, ~final_state, ~decision_hour,
    "A", 1L, "ruled_out", 2L,
    "B", 0L, "ruled_out", 3L,
    "C", 1L, "not_eligible", 1L
  )
  cohort <- tibble::tibble(
    study_id = c("A", "A", "B", "C"),
    picu_adm_date_time = as.POSIXct(
      c("2026-01-01 08:00:00", "2026-01-01 08:00:00",
        "2026-01-02 09:00:00", "2026-01-03 10:00:00"),
      tz = "America/Denver"
    )
  )

  result <- create_false_negative_chart_review_table(
    decisions,
    cohort,
    expected_n = 1L
  )

  expect_named(result, c(
    "study_id", "picu_adm_date_time", "hours_after_picu_admission"
  ))
  expect_equal(result$study_id, "A")
  expect_equal(result$hours_after_picu_admission, 2)
  expect_equal(result$picu_adm_date_time, cohort$picu_adm_date_time[[1]])
})

test_that("chart-review table detects unsafe or unexpected rosters", {
  decisions <- tibble::tibble(
    study_id = "A", sbi_present = 1L,
    final_state = "ruled_out", decision_hour = 2L
  )
  duplicate_cohort <- tibble::tibble(
    study_id = c("A", "A"),
    picu_adm_date_time = as.POSIXct(
      c("2026-01-01 08:00:00", "2026-01-01 09:00:00"), tz = "UTC"
    )
  )

  expect_error(
    create_false_negative_chart_review_table(decisions, duplicate_cohort),
    "Multiple PICU admission times"
  )
  expect_error(
    create_false_negative_chart_review_table(
      decisions,
      duplicate_cohort[1, ],
      expected_n = 18L
    ),
    "Expected 18"
  )
})
