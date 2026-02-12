library(testthat)

validate_predictor_schema <- function(df, required_predictors, dataset_label) {
  missing_predictors <- setdiff(required_predictors, colnames(df))
  if (length(missing_predictors) > 0) {
    stop(
      paste0(
        dataset_label,
        " is missing required predictors: ",
        paste(missing_predictors, collapse = ", ")
      )
    )
  }

  extra_predictors <- setdiff(colnames(df), required_predictors)
  list(missing = missing_predictors, extra = extra_predictors)
}

test_that("validate_predictor_schema passes when all required predictors are present", {
  df <- data.frame(a = 1, b = 2, c = 3)
  result <- validate_predictor_schema(df, c("a", "b"), "toy_df")

  expect_length(result$missing, 0)
  expect_true(all(c("c") %in% result$extra))
})

test_that("validate_predictor_schema fails when required predictors are missing", {
  df <- data.frame(a = 1, c = 3)

  expect_error(
    validate_predictor_schema(df, c("a", "b"), "toy_df"),
    "toy_df is missing required predictors: b"
  )
})
