source("plot_save_helpers.R")
### Evaluate how predictor count affects AUROC for SBI RF model ###

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
  library(doParallel)
  library(pROC)
  library(fst)
  library(tidyverse)
  library(tidyr)
  library(readr)
  library(dplyr)
  library(readxl)
  library(purrr)
  library(chron)
  library(stringr)
  library(stringi)
  library(stats)
  library(tibble)
  library(robustbase) #For colMedians function
  library(DescTools)
  library(rriskDistributions)
  library(tigerstats)
  library(fst)
  library(rpart)
  library(ranger)
  library(caTools)
  library(caret)
  library(pROC)
  library(pccc)
  library(stratifyR)
  library(glmnet)
  library(regclass)
  library(devtools)
  library(recipes)
  library(ggpubr)
  library(e1071)
  library(mltools)
  library(data.table)
  library(kernlab)
  library(MLmetrics)
  library(OptimalCutpoints)
  library(qwraps2)
})

# Load in dataframe
temp_model_data_df <- read.fst(path = "~/sbi_blake/jan_25_23_model_data_df.fst")
model_data_df <- temp_model_data_df

# Filter out post-Beaker patients
t_beaker <- as.POSIXct("2019-10-31 00:00:00 UTC", tz = "UTC")
vps_pt_list_11_to_17 <- read_excel("/phi/sbi/sbi_data/vps_full_admit_list.xlsx")
vps_pt_list_18_to_20 <- read_excel(vps_file_path_10_yr, sheet = "Pt List")
vps_pt_list_full <- bind_rows(vps_pt_list_11_to_17, vps_pt_list_18_to_20)
vps_pt_list_full <- vps_pt_list_full %>% mutate(mrn = as.character(vps_pt_list_full$mrn),
                                                case_id = as.character(vps_pt_list_full$case_id),
                                                har = as.character(vps_pt_list_full$har),
                                                picu_adm_date_time = as.POSIXct(vps_pt_list_full$picu_adm_date_time, tz = "GMT"),
                                                picu_dc_date_time = as.POSIXct(vps_pt_list_full$picu_dc_date_time, tz = "GMT"))

case_id_and_picu_adm_date_time <- vps_pt_list_full %>% dplyr::select(case_id, picu_adm_date_time) %>% distinct()

# Create base dataframe for all pre-Beaker PICU encounters
rf_df_pre_beaker <- model_data_df %>%
  left_join(case_id_and_picu_adm_date_time, by = "case_id") %>%
  filter(picu_adm_date_time < t_beaker) %>%
  dplyr::select(-all_of(c("picu_adm_date_time")))

cat(sprintf("Cohort size after pre-Beaker filter (all abx exposure statuses): %s\n", format(nrow(rf_df_pre_beaker), big.mark = ",")))

run_rf_predictor_sweep <- function(rf_df, cohort_label, workers, positive_class = "no") {
  bad_cols <- c("sbi_pneumonia", "sbi", "blood", "abd", "cns", "genital", "nos", "resp", "stool", "tissue", "urine", "sbi_cx_neg_sepsis", "virus_24", "abx_exp", "death_date")
  rf_df <- rf_df %>% dplyr::select(-all_of(bad_cols))

  # Fix classes
  rf_df$sbi_present <- as.factor(rf_df$sbi_present)

  # Remove impossible values
  rf_df <- rf_df %>% filter(max_sbp <= 100)

  cat(sprintf("\n=== Running cohort: %s (n=%s after cleaning) ===\n", cohort_label, format(nrow(rf_df), big.mark = ",")))

  # 0) Train / test split
  set.seed(32313)
  train_size <- floor(nrow(rf_df) * 0.8)
  idx <- sample.int(n = nrow(rf_df), size = train_size, replace = FALSE)
  train_df <- rf_df[idx, ]
  test_df <- rf_df[-idx, ]

  train_df$sbi_present <- relevel(train_df$sbi_present, ref = positive_class)
  test_df$sbi_present <- factor(test_df$sbi_present, levels = levels(train_df$sbi_present))

  # Never allow encounter identifier to be used as a predictor
  train_df <- train_df[, setdiff(names(train_df), "case_id"), drop = FALSE]
  test_df <- test_df[, setdiff(names(test_df), "case_id"), drop = FALSE]

  # 1) Repeated CV for tuning on train only
  ctrl_repeated <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 2,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    allowParallel = TRUE,
    verboseIter = FALSE
  )

  # 2) Tune RF
  set.seed(32313)
  rf_tuned <- caret::train(
    sbi_present ~ .,
    data = train_df,
    method = "ranger",
    metric = "ROC",
    trControl = ctrl_repeated,
    tuneLength = 20,
    importance = "none",
    respect.unordered.factors = "order",
    num.threads = 1,
    na.action = na.omit
  )

  best_params <- rf_tuned$bestTune
  cat(sprintf("Best tuned parameters for %s:\n", cohort_label))
  print(best_params)

  # 3) Permutation importance ranking on full train
  set.seed(32313)
  rf_importance <- ranger(
    formula = sbi_present ~ ., data = train_df,
    mtry = best_params$mtry,
    splitrule = as.character(best_params$splitrule),
    min.node.size = best_params$min.node.size,
    num.trees = 500, probability = TRUE,
    importance = "permutation", respect.unordered.factors = "order",
    num.threads = workers, na.action = "na.omit"
  )

  importance_df <- data.frame(
    feature = names(rf_importance$variable.importance),
    importance = as.numeric(rf_importance$variable.importance),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(importance))

  predictors_ranked <- importance_df$feature
  p <- length(predictors_ranked)

  # 4) For k = 1..p, fit RF with top-k features via repeated CV
  roc_by_k <- numeric(p)

  for (k in seq_len(p)) {
    top_k <- predictors_ranked[seq_len(k)]
    train_k <- train_df[, c("sbi_present", top_k), drop = FALSE]

    tune_k <- best_params
    tune_k$mtry <- min(tune_k$mtry, k)

    set.seed(32313 + k)
    fit_k <- train(
      sbi_present ~ ., data = train_k, method = "ranger", metric = "ROC",
      trControl = ctrl_repeated, tuneGrid = tune_k,
      importance = "none", respect.unordered.factors = "order",
      num.threads = 1, na.action = na.omit
    )

    roc_by_k[k] <- max(fit_k$results$ROC, na.rm = TRUE)
    cat(sprintf("[%s] k=%d/%d, CV AUROC=%.4f\n", cohort_label, k, p, roc_by_k[k]))
  }

  performance_by_k <- data.frame(k = seq_len(p), cv_auroc = roc_by_k)

  max_auc <- max(performance_by_k$cv_auroc, na.rm = TRUE)
  k_max_auc <- performance_by_k$k[which.max(performance_by_k$cv_auroc)]

  # 5) Plot

  # Create cleaner label for plot title
  plot_cohort_label <- dplyr::case_when(
    cohort_label == "abx_exp == 0" ~ "antibiotic unexposed",
    cohort_label == "abx_exp == 1" ~ "antibiotic exposed",
    TRUE ~ cohort_label
  )

  p_cv_auroc_by_predictor_count <- ggplot(performance_by_k, aes(x = k, y = cv_auroc)) +
      geom_line(linewidth = 0.8, color = "red") +
      geom_point(size = 2, color = "red") +
      geom_vline(
        xintercept = k_max_auc,
        linetype = "dashed",
        color = "steelblue",
        linewidth = 0.8
      ) +
      geom_hline(
        yintercept = max_auc,
        linetype = "dotted",
        color = "darkgray",
        linewidth = 0.8
      ) +
      scale_x_continuous(
        breaks = seq(
          from = 0,
          to = max(performance_by_k$k, na.rm = TRUE),
          by = 10
        )
      ) +
      labs(
        title = paste0(
          "CV AUROC vs Number of Predictors (",
          plot_cohort_label,
          ")"
        ),
        x = "Number of predictors (top-k by permutation importance)",
        y = "Cross-validated AUROC"
      ) +
      coord_cartesian(
        ylim = range(performance_by_k$cv_auroc, na.rm = TRUE)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      )

  print(p_cv_auroc_by_predictor_count)
  save_aim1_plot(
    p_cv_auroc_by_predictor_count,
    paste0("cv_auroc_by_predictor_count_", gsub("[^A-Za-z0-9]+", "_", tolower(plot_cohort_label)), ".tiff")
  )

  # 6) Parsimony rule: smallest k within 99% of max AUROC
  threshold_auc <- 0.99 * max_auc
  selected_k <- min(performance_by_k$k[performance_by_k$cv_auroc >= threshold_auc])
  selected_features <- predictors_ranked[seq_len(selected_k)]

  cat(sprintf("[%s] Max CV AUROC: %.4f\n", cohort_label, max_auc))
  cat(sprintf("[%s] Parsimony threshold (99%% of max): %.4f\n", cohort_label, threshold_auc))
  cat(sprintf("[%s] Selected k: %d\n", cohort_label, selected_k))

  auroc_by_predictor_table <- performance_by_k %>%
    dplyr::transmute(n_predictors = k, cv_auroc = round(cv_auroc, 4))

  View(auroc_by_predictor_table)

  # 7) Refit final model on full train with selected k
  final_train <- train_df[, c("sbi_present", selected_features), drop = FALSE]
  final_tune <- best_params
  final_tune$mtry <- min(final_tune$mtry, selected_k)

  set.seed(32313)
  final_rf <- train(
    sbi_present ~ ., data = final_train, method = "ranger", metric = "ROC",
    trControl = trainControl(method = "none", classProbs = TRUE),
    tuneGrid = final_tune,
    importance = "none", respect.unordered.factors = "order",
    num.threads = workers, na.action = na.omit
  )

  # 8) Evaluate once on test_df (AUROC + CI)
  final_test <- test_df[, c("sbi_present", selected_features), drop = FALSE]
  probs <- predict(final_rf, newdata = final_test, type = "prob")[[positive_class]]

  roc_obj <- pROC::roc(
    response = final_test$sbi_present,
    predictor = probs,
    levels = rev(levels(final_test$sbi_present)),
    direction = "<",
    ci = TRUE
  )

  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))

  cat(sprintf("[%s] Final test AUROC: %.4f\n", cohort_label, auc_val))
  cat(sprintf("[%s] 95%% CI: [%.4f, %.4f]\n", cohort_label, auc_ci[1], auc_ci[3]))

  list(
    cohort = cohort_label,
    tuned_model = rf_tuned,
    importance = importance_df,
    performance_by_k = performance_by_k,
    selected_k = selected_k,
    selected_features = selected_features,
    final_model = final_rf,
    test_auc = auc_val,
    test_auc_ci = auc_ci
  )
}

#-----------------------------
# Parallel backend (n - 1 cores)
#-----------------------------
n_cores <- parallel::detectCores(logical = TRUE)
workers <- max(1, n_cores - 1)
cl <- parallel::makePSOCKcluster(workers)
doParallel::registerDoParallel(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}, add = TRUE)

# Run separate analyses by antibiotic exposure status
rf_df_abx_unexposed <- rf_df_pre_beaker %>% filter(abx_exp == "0")
rf_df_abx_exposed <- rf_df_pre_beaker %>% filter(abx_exp == "1")

results_abx_unexposed <- run_rf_predictor_sweep(
  rf_df = rf_df_abx_unexposed,
  cohort_label = "abx_exp == 0",
  workers = workers
)

results_abx_exposed <- run_rf_predictor_sweep(
  rf_df = rf_df_abx_exposed,
  cohort_label = "abx_exp == 1",
  workers = workers
)

best_params <- rf_tuned$bestTune
cat("Best tuned parameters:\n")
print(best_params)

# Optional objects for downstream use
results <- list(
  abx_unexposed = results_abx_unexposed,
  abx_exposed = results_abx_exposed
)
