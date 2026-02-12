# Script to run retrospective models on prospective data
# and evaluate performance restricted to prospective patients with suspected infection.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pROC)
  library(yardstick)
  library(gtsummary)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1) Characteristics: prospective suspected infection vs no suspected infection
# -----------------------------------------------------------------------------

pros_characteristics <- bind_rows(
  pros_no_abx_1st_infxn %>% mutate(abx_exposure = "No antibiotics"),
  pros_yes_abx_1st_infxn %>% mutate(abx_exposure = "Antibiotics")
) %>%
  mutate(
    suspected_infection = factor(
      if_else(as.integer(suspected_infection) == 1L, "Suspected infection", "No suspected infection"),
      levels = c("Suspected infection", "No suspected infection")
    ),
    is_female = as.character(is_female),
    is_female = case_when(
      is_female %in% c("1", "TRUE", "T") ~ "Female",
      is_female %in% c("0", "FALSE", "F") ~ "Male",
      TRUE ~ NA_character_
    ),
    is_female = factor(is_female)
  )

characteristic_vars <- c(
  "abx_exposure", "age", "is_female", "race", "ethnicity", "pccc",
  "malignancy_pccc", "los_before_icu_days", "imv_at_picu_adm", "sbi_present"
)
characteristic_vars <- intersect(characteristic_vars, colnames(pros_characteristics))

table1_pros_suspicion <- pros_characteristics %>%
  select(all_of(c("suspected_infection", characteristic_vars))) %>%
  tbl_summary(
    by = suspected_infection,
    statistic = list(
      all_continuous() ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing = "ifany"
  ) %>%
  add_p() %>%
  modify_header(label ~ "**Characteristic**") %>%
  bold_labels()

print(table1_pros_suspicion)

# -----------------------------------------------------------------------------
# 2) Restrict prospective cohorts to suspected infection only
# -----------------------------------------------------------------------------
pros_no_abx_susp <- pros_no_abx_1st_infxn %>% filter(as.integer(suspected_infection) == 1L)
pros_yes_abx_susp <- pros_yes_abx_1st_infxn %>% filter(as.integer(suspected_infection) == 1L)

# For yes-abx model, align names to expected predictors
predictors_abx[predictors_abx == "hematocrit_mean"] <- "hematocrit_blood"
pros_yes_temp <- pros_yes_abx_susp %>%
  rename(hematocrit_blood = hematocrit_mean, last_fio2 = fio2_last, max_dbp = dbp_max)

# Score suspected-infection subsets with retrospective RF models
rf_pred_prob_pros <- predict(rf_model, pros_no_abx_susp[, all_of(predictors)], type = "prob")
rf_pred_prob_pros_abx <- predict(rf_model_abx, pros_yes_temp[, all_of(predictors_abx)], type = "prob")

pros_no_eval <- tibble(
  cohort = "Prospective suspected infection • Abx-",
  truth_num = as.integer(pros_no_abx_susp$sbi_present),
  truth = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes", "no")),
  p_yes = rf_pred_prob_pros[, "yes"]
)

pros_yes_eval <- tibble(
  cohort = "Prospective suspected infection • Abx+",
  truth_num = as.integer(pros_yes_temp$sbi_present),
  truth = factor(if_else(truth_num == 1L, "yes", "no"), levels = c("yes", "no")),
  p_yes = rf_pred_prob_pros_abx[, "yes"]
)

# -----------------------------------------------------------------------------
# 3) Metrics helper (AUROC, AUPRC, NPV, false negatives)
# -----------------------------------------------------------------------------

metric_summary <- function(eval_df, threshold) {
  pred_neg <- eval_df$p_yes < threshold
  actual_pos <- eval_df$truth_num == 1L
  actual_neg <- eval_df$truth_num == 0L

  tn <- sum(pred_neg & actual_neg, na.rm = TRUE)
  fn <- sum(pred_neg & actual_pos, na.rm = TRUE)

  tibble(
    n = nrow(eval_df),
    prevalence = mean(actual_pos, na.rm = TRUE),
    auroc = as.numeric(pROC::auc(pROC::roc(eval_df$truth_num, eval_df$p_yes, quiet = TRUE, direction = "<"))),
    auprc = yardstick::pr_auc(eval_df, truth = truth, p_yes) %>% pull(.estimate),
    threshold = threshold,
    npv = if_else((tn + fn) > 0, tn / (tn + fn), as.numeric(NA)),
    false_negative_n = fn,
    false_negative_rate_among_sbi = if_else(sum(actual_pos, na.rm = TRUE) > 0, fn / sum(actual_pos, na.rm = TRUE), as.numeric(NA))
  )
}

pros_metrics_no_abx <- metric_summary(pros_no_eval, threshold = 0.05) %>%
  mutate(cohort = "Prospective suspected infection • Abx-")

pros_metrics_yes_abx <- metric_summary(pros_yes_eval, threshold = 0.074) %>%
  mutate(cohort = "Prospective suspected infection • Abx+")

pros_metrics_suspected <- bind_rows(pros_metrics_no_abx, pros_metrics_yes_abx) %>%
  select(cohort, everything())

print(pros_metrics_suspected)

# -----------------------------------------------------------------------------
# 4) Build retrospective comparison sets (test set performance)
# -----------------------------------------------------------------------------
retro_no_eval <- test_pr_df %>%
  transmute(
    cohort = "Retrospective test set • Abx-",
    truth = factor(as.character(truth), levels = c("yes", "no")),
    truth_num = as.integer(truth == "yes"),
    p_yes = as.numeric(p_yes)
  )

retro_yes_eval <- test_pr_df_abx %>%
  transmute(
    cohort = "Retrospective test set • Abx+",
    truth = factor(as.character(truth), levels = c("yes", "no")),
    truth_num = as.integer(truth == "yes"),
    p_yes = as.numeric(p_yes)
  )

# -----------------------------------------------------------------------------
# 5) ROC + PR curve utilities and plots (retro test vs prospective suspected)
# -----------------------------------------------------------------------------
roc_df <- function(df, cohort_name) {
  roc_obj <- pROC::roc(df$truth_num, df$p_yes, quiet = TRUE, direction = "<")
  tibble(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    cohort = cohort_name
  )
}

pr_df <- function(df, cohort_name) {
  out <- yardstick::pr_curve(df, truth = truth, p_yes)
  out %>% transmute(recall = recall, precision = precision, cohort = cohort_name)
}

roc_no <- bind_rows(
  roc_df(retro_no_eval, "Retrospective test set"),
  roc_df(pros_no_eval, "Prospective suspected infection")
)

roc_yes <- bind_rows(
  roc_df(retro_yes_eval, "Retrospective test set"),
  roc_df(pros_yes_eval, "Prospective suspected infection")
)

pr_no <- bind_rows(
  pr_df(retro_no_eval, "Retrospective test set"),
  pr_df(pros_no_eval, "Prospective suspected infection")
)

pr_yes <- bind_rows(
  pr_df(retro_yes_eval, "Retrospective test set"),
  pr_df(pros_yes_eval, "Prospective suspected infection")
)

plot_roc_no <- ggplot(roc_no, aes(x = fpr, y = tpr, color = cohort)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  coord_equal() +
  labs(
    title = "ROC: Abx-unexposed model",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_bw()

plot_roc_yes <- ggplot(roc_yes, aes(x = fpr, y = tpr, color = cohort)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  coord_equal() +
  labs(
    title = "ROC: Abx-exposed model",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_bw()

plot_pr_no <- ggplot(pr_no, aes(x = recall, y = precision, color = cohort)) +
  geom_line(linewidth = 1) +
  labs(
    title = "PR: Abx-unexposed model",
    x = "Recall",
    y = "Precision",
    color = NULL
  ) +
  theme_bw()

plot_pr_yes <- ggplot(pr_yes, aes(x = recall, y = precision, color = cohort)) +
  geom_line(linewidth = 1) +
  labs(
    title = "PR: Abx-exposed model",
    x = "Recall",
    y = "Precision",
    color = NULL
  ) +
  theme_bw()

rf_curve_comparison_plot <- (plot_roc_no | plot_roc_yes) / (plot_pr_no | plot_pr_yes)
print(rf_curve_comparison_plot)
