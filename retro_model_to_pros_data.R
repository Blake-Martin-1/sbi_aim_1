# Script to run retrospective model on prospective input data

# Now test the no_abx model on prospective 2hr AU cohort
rf_pred_prob_pros <- predict(rf_model, pros_no_abx_1st_infxn[, all_of(c(predictors))], type = "prob")
pros_auroc_no_abx <- colAUC(rf_pred_prob_pros, pros_no_abx_1st_infxn$sbi_present, plotROC = FALSE) # AUC of 0.736 in prospective test set

# Now test the yes_abx model on prospective 2hr AU cohort
# First fix predictors_abx to match names in pros df
predictors_abx[predictors_abx == "hematocrit_mean"] <- "hematocrit_blood"

pros_yes_temp <- pros_yes_abx_1st_infxn %>% rename(hematocrit_blood = hematocrit_mean, last_fio2 = fio2_last, max_dbp = dbp_max)

rf_pred_prob_pros_abx <- predict(rf_model_abx, (pros_yes_temp[, all_of(c(predictors_abx))]), type = "prob")
pros_auroc_abx_yes_abx <- colAUC(rf_pred_prob_pros_abx, pros_yes_temp$sbi_present, plotROC = FALSE) # AUC of 0.750 in prospective test set



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
  truth = factor(ifelse(pros_no_abx_1st_infxn$sbi_present == 1, "yes", "no"),
                 levels = c("yes","no")),
  p_yes = rf_pred_prob_pros[, "yes"]
)

pros_auprc_no_abx <- pr_auc(pros_noabx_pr_df, truth = truth, p_yes) %>%
  pull(.estimate)

print(paste0("Prospective no-abx AUPRC = ", signif(pros_auprc_no_abx, 4)))


# ---- Yes-abx prospective AUPRC ----
stopifnot(all(c("yes","no") %in% colnames(rf_pred_prob_pros_abx)))
stopifnot(all(pros_yes_temp$sbi_present %in% c(0, 1)))

pros_yesabx_pr_df <- tibble(
  truth = factor(ifelse(pros_yes_temp$sbi_present == 1, "yes", "no"),
                 levels = c("yes","no")),
  p_yes = rf_pred_prob_pros_abx[, "yes"]
)

pros_auprc_yes_abx <- pr_auc(pros_yesabx_pr_df, truth = truth, p_yes) %>%
  pull(.estimate)

print(paste0("Prospective yes-abx AUPRC = ", signif(pros_auprc_yes_abx, 4)))
