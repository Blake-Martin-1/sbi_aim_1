library(testthat)

source("../../prospective_sbi_type_breakdown.R")

test_that("SBI types report n (%) using each exposure stratum denominator", {
  origin <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  unexposed <- tibble::tibble(
    study_id = c("u1", "u2", "u3"),
    mrn = c("1", "2", "3"),
    picu_adm_date_time = origin,
    ever_cx_neg_sepsis = c(0, 1, 0),
    pna_1_0 = c(1, 0, 0)
  )
  exposed <- tibble::tibble(
    study_id = c("e1", "u1", "e2"),
    mrn = c("4", "1", "5"),
    picu_adm_date_time = origin,
    ever_cx_neg_sepsis = 0,
    pna_1_0 = 0
  )
  micro <- tibble::tibble(
    mrn = c("1", "1", "1", "4", "5"),
    specimen_source = c(
      "Blood, Venous", "Blood, Venous", "Urine, Clean Catch",
      "Nasopharyngeal\u00a0Swab", "Blood, Venous"
    ),
    time_obtained = origin + c(0, 60, 120, 0, 25 * 60 * 60)
  )

  result <- summarize_prospective_sbi_types(micro, unexposed, exposed)
  summary <- result$sbi_type_summary

  expect_equal(
    summary$antibiotic_unexposed[summary$broad_category == "blood"],
    "1 (33.3%)"
  )
  expect_equal(
    summary$antibiotic_exposed[summary$broad_category == "blood"],
    "1 (33.3%)"
  )
  expect_equal(
    summary$antibiotic_exposed[summary$broad_category == "respiratory"],
    "1 (33.3%)"
  )
  expect_equal(
    summary$antibiotic_unexposed[summary$broad_category == "culture_negative_sepsis"],
    "1 (33.3%)"
  )
  expect_equal(
    summary$antibiotic_unexposed[summary$broad_category == "bacterial_pneumonia"],
    "1 (33.3%)"
  )
  expect_equal(names(summary), c(
    "broad_category", "antibiotic_unexposed", "antibiotic_exposed"
  ))

  # u1 has blood, urine, and VPS bacterial pneumonia and must remain in all
  # three categories rather than being forced into one mutually exclusive type.
  expect_setequal(
    result$encounter_sbi_types$broad_category[result$encounter_sbi_types$study_id == "u1"],
    c("blood", "gu", "bacterial_pneumonia")
  )
})

test_that("encounter audit rows retain distinct sources and broad types", {
  origin <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  unexposed <- tibble::tibble(
    study_id = "u1", mrn = "1", picu_adm_date_time = origin,
    ever_cx_neg_sepsis = 0, pna_1_0 = 0
  )
  exposed <- tibble::tibble(
    study_id = character(), mrn = character(),
    picu_adm_date_time = as.POSIXct(character()),
    ever_cx_neg_sepsis = numeric(), pna_1_0 = numeric()
  )
  micro <- tibble::tibble(
    mrn = c("1", "1"),
    specimen_source = c("Craniotomy", "new source"),
    time_obtained = origin
  )

  details <- summarize_prospective_sbi_types(micro, unexposed, exposed)$encounter_sbi_types

  expect_setequal(details$broad_category, c("cns", "unknown"))
  expect_equal(unique(details$study_id), "u1")
})

test_that("snake_case sources from pros_micro_slim map to anatomical categories", {
  origin <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  unexposed <- tibble::tibble(
    study_id = "u1", mrn = "1", picu_adm_date_time = origin,
    ever_cx_neg_sepsis = 0, pna_1_0 = 0
  )
  exposed <- unexposed[0, ]
  micro <- tibble::tibble(
    mrn = "1",
    specimen_source = c(
      "peripheral_venipuncture", "peripheral_piv_start", "catheter_port",
      "catheter_picc", "urine_clean_catch", "lumbar_puncture",
      "nasopharyngeal_swab"
    ),
    time_obtained = origin
  )

  details <- summarize_prospective_sbi_types(
    micro, unexposed, exposed
  )$encounter_sbi_types

  expect_setequal(
    details$broad_category,
    c("blood", "gu", "cns", "respiratory")
  )
  expect_false(any(details$broad_category == "unknown"))
  expect_true("peripheral_venipuncture" %in% details$specimen_source)
})

test_that("a patient in both exposure strata is counted once in each", {
  origin <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  cohort <- tibble::tibble(
    study_id = "same", mrn = "1", picu_adm_date_time = origin,
    ever_cx_neg_sepsis = 0, pna_1_0 = 0
  )
  micro <- tibble::tibble(
    mrn = "1", specimen_source = "Blood, Venous", time_obtained = origin
  )

  summary <- summarize_prospective_sbi_types(micro, cohort, cohort)$sbi_type_summary

  blood <- summary[summary$broad_category == "blood", ]
  expect_equal(blood$antibiotic_unexposed, "1 (100.0%)")
  expect_equal(blood$antibiotic_exposed, "1 (100.0%)")
})
