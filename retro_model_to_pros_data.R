# Script to run retrospective model on prospective input data

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

bootstrap_metric_ci <- function(truth_case, score_case, metric = c("auroc", "auprc"), n_boot = 1000) {
  metric <- match.arg(metric)
  eval_df <- tibble::tibble(truth_case = truth_case, score_case = score_case)

  if (nrow(eval_df) == 0 || dplyr::n_distinct(eval_df$truth_case) < 2) {
    return(c(NA_real_, NA_real_))
  }

  boot_vals <- replicate(n_boot, {
    idx <- sample.int(nrow(eval_df), size = nrow(eval_df), replace = TRUE)
    b <- eval_df[idx, , drop = FALSE]
    if (dplyr::n_distinct(b$truth_case) < 2) return(NA_real_)

    if (metric == "auroc") {
      as.numeric(pROC::auc(pROC::roc(b$truth_case, b$score_case, quiet = TRUE, direction = "<")))
    } else {
      PRROC::pr.curve(
        scores.class0 = b$score_case[b$truth_case == 1],
        scores.class1 = b$score_case[b$truth_case == 0],
        curve = FALSE
      )$auc.integral
    }
  })

  boot_vals <- boot_vals[is.finite(boot_vals)]
  if (length(boot_vals) < 20) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE))
}

fmt_ci <- function(est, ci) sprintf("%.2f (%.2f - %.2f)", round(est, 2), round(ci[1], 2), round(ci[2], 2))

# Now test the no_abx model on prospective 2hr AU cohort
validate_predictor_schema(pros_no_abx_1st_infxn, predictors, "pros_no_abx_1st_infxn")
rf_pred_prob_pros <- predict(rf_model, pros_no_abx_1st_infxn[, all_of(c(predictors))], type = "prob")
pros_auroc_no_abx <- colAUC(1 - rf_pred_prob_pros, 1 - pros_no_abx_1st_infxn$sbi_present, plotROC = FALSE) # AUC using SBI-negative as cases
pros_auroc_no_abx


# Create dataframe to view probabilities and outcomes
pred_df_pros_no_abx <- pros_no_abx_1st_infxn %>%
  transmute(
    rowIndex = row_number(),
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros[, "yes"]
  )

pros_auroc_no_abx_est <- as.numeric(
  pROC::auc(
    pROC::roc(
      response = as.integer(pred_df_pros_no_abx$sbi_present == 0L),
      predictor = 1 - pred_df_pros_no_abx$pred_prob_yes,
      quiet = TRUE,
      direction = "<"
    )
  )
)
pros_auroc_no_abx_ci <- bootstrap_metric_ci(
  truth_case = as.integer(pred_df_pros_no_abx$sbi_present == 0L),
  score_case = 1 - pred_df_pros_no_abx$pred_prob_yes,
  metric = "auroc",
  n_boot = 1000
)
print(paste0("Prospective no-abx AUROC (95% CI) = ", fmt_ci(pros_auroc_no_abx_est, pros_auroc_no_abx_ci)))




# Now test the yes_abx model on prospective 2hr AU cohort
# First fix predictors_abx to match names in pros df
predictors_abx[predictors_abx == "hematocrit_mean"] <- "hematocrit_blood"

pros_yes_temp <- pros_yes_abx_1st_infxn %>% rename(hematocrit_blood = hematocrit_mean, last_fio2 = fio2_last, max_dbp = dbp_max)

validate_predictor_schema(pros_yes_temp, predictors_abx, "pros_yes_temp")
rf_pred_prob_pros_abx <- predict(rf_model_abx, (pros_yes_temp[, all_of(c(predictors_abx))]), type = "prob")
pros_auroc_abx_yes_abx <- colAUC(1 - rf_pred_prob_pros_abx, 1 - pros_yes_temp$sbi_present, plotROC = FALSE) # AUC using SBI-negative as cases
pros_auroc_abx_yes_abx


# Now view the predictions along with probabilities if you want
pred_df_pros_abx <- pros_yes_temp %>%
  transmute(
    sbi_present = sbi_present,
    pred_prob_yes = rf_pred_prob_pros_abx[, "yes"]
  )

# Get AUROC
roc_obj_pros_yes_abx <- roc(
  response  = pred_df_pros_abx$sbi_present,
  predictor = 1 - pred_df_pros_abx$pred_prob_yes,
  levels    = c(1, 0),   # 0 (SBI-negative) is the positive class
  direction = "<"
)

pros_auroc_yes_abx <- as.numeric(auc(roc_obj_pros_yes_abx))
pros_auroc_yes_abx_ci <- bootstrap_metric_ci(
  truth_case = as.integer(pred_df_pros_abx$sbi_present == 0L),
  score_case = 1 - pred_df_pros_abx$pred_prob_yes,
  metric = "auroc",
  n_boot = 1000
)
print(paste0("Prospective yes-abx AUROC (95% CI) = ", fmt_ci(pros_auroc_yes_abx, pros_auroc_yes_abx_ci)))


# =========================
# Prospective AUPRC (no_abx + yes_abx)
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(yardstick)
})

# ---- No-abx prospective AUPRC ----
stopifnot(all(c("yes","no") %in% colnames(rf_pred_prob_pros)))   # caret probs
stopifnot(all(pros_no_abx_1st_infxn$sbi_present %in% c(0, 1)))   # numeric truth

pros_noabx_pr_df <- tibble(
  truth = factor(ifelse(pros_no_abx_1st_infxn$sbi_present == 0, "yes", "no"),
                 levels = c("yes","no")),
  p_yes = 1 - rf_pred_prob_pros[, "yes"]
)

pros_auprc_no_abx <- pr_auc(pros_noabx_pr_df, truth = truth, p_yes) %>%
  pull(.estimate)
pros_auprc_no_abx_ci <- bootstrap_metric_ci(
  truth_case = as.integer(pros_no_abx_1st_infxn$sbi_present == 0L),
  score_case = 1 - rf_pred_prob_pros[, "yes"],
  metric = "auprc",
  n_boot = 1000
)

print(paste0("Prospective no-abx AUPRC (95% CI) = ", fmt_ci(pros_auprc_no_abx, pros_auprc_no_abx_ci)))


# ---- Yes-abx prospective AUPRC ----
stopifnot(all(c("yes","no") %in% colnames(rf_pred_prob_pros_abx)))
stopifnot(all(pros_yes_temp$sbi_present %in% c(0, 1)))

pros_yesabx_pr_df <- tibble(
  truth = factor(ifelse(pros_yes_temp$sbi_present == 0, "yes", "no"),
                 levels = c("yes","no")),
  p_yes = 1 - rf_pred_prob_pros_abx[, "yes"]
)

pros_auprc_yes_abx <- pr_auc(pros_yesabx_pr_df, truth = truth, p_yes) %>%
  pull(.estimate)
pros_auprc_yes_abx_ci <- bootstrap_metric_ci(
  truth_case = as.integer(pros_yes_temp$sbi_present == 0L),
  score_case = 1 - rf_pred_prob_pros_abx[, "yes"],
  metric = "auprc",
  n_boot = 1000
)

print(paste0("Prospective yes-abx AUPRC (95% CI) = ", fmt_ci(pros_auprc_yes_abx, pros_auprc_yes_abx_ci)))
