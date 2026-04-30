### Code to create new model from prospective data ###


pros_one_model <- read_csv(file = "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_10_14_25.csv")
# pros_all <- read_csv(file = "/phi/sbi/sbi_blake/pros_all_just_b4_modeling_1_15_26_all_models.csv")

### Train new models using prospective data ###
# slim down to just the study_id, predictor info, and sbi outcome for model training


model_slim <- pros_one_model %>% dplyr::select(-pat_enc_csn_id, -picu_adm_date_time, -intime, -outtime,
                                               -contact_creation_to_admission_delta, -height, -height_date, -origin, -valid_start_instant, -valid_end_instant, -cx_neg_sepsis,
                                               -pna_1_0, -bcx_sent, -micro_sbi_1_0, -cnss_true, -ever_cx_neg_sepsis, -hsp_account_id, -pat_mrn_id, -race, -ethnicity,
                                               -language, -insurance_type, -icu_start_instant, -old_race, -prior_unit, -weight, -weight_date, -is_interventional_radiology,
                                               -abx_duration_after_score
)

raw_names <- colnames(model_slim)[str_ends(colnames(model_slim), "_raw")]

model_slim <- model_slim %>% dplyr::select(-all_of(raw_names))

# Re-order
model_slim <- model_slim %>% relocate(study_id, sbi_present, everything())

# Remove duplicate rows
model_slim <- model_slim %>% distinct() #62,814


# Now some QC to see IDs where SBI-present changes:
# IDs with mixed SBI labels
mixed_ids <- model_slim %>%
  group_by(study_id) %>%
  filter(n_distinct(sbi_present) > 1) %>%      # TRUE if both 0 and 1 seen
  ungroup() # currently no rows

## Remove bad columns
predictors <- setdiff(names(model_slim), c("study_id", "sbi_present", "SBI", "abx_exp", "rowid", "hours_since_picu_adm"))

nzv <- nearZeroVar(model_slim[ , predictors])
if(length(nzv) > 0) {
  predictors <- predictors[-nzv]
} # gets rid of 28 predictors

bad_cols <- names(model_slim)[
  sapply(model_slim[ , predictors], function(x) any(is.na(x) | is.infinite(x)))
] # no bad columns

# At this point there are 69 predictors

#Set seed for reproducibility
set.seed(2025)

# Collapse to one row per admission to obtain PICU encounter-level labels
admission_df <- model_slim %>%
  group_by(study_id) %>% slice(1) %>% ungroup() %>%
  mutate(SBI_enc = ifelse(sbi_present == 1, "pos", "neg"))

# Shuffle IDs
all_ids <- unique(admission_df$study_id)
shuffled_ids <- sample(all_ids)

# 60/20/20 split fir train, validation, and test sets
n_total <- length(shuffled_ids)
n_train <- floor(0.60 * n_total)
n_valid <- floor(0.20 * n_total)

train_ids <- shuffled_ids[1:n_train]
valid_ids <- shuffled_ids[(n_train + 1):(n_train + n_valid)]
test_ids  <- shuffled_ids[(n_train + n_valid + 1):n_total]

# Save splits to enable ability to rebuild model in Epic python environment
splits_list <- list(
  seed = 2025L,
  train_ids = as.character(train_ids),
  valid_ids = as.character(valid_ids),
  test_ids  = as.character(test_ids)
)

# Build dataframes with all the predictors and SBI outcome label using the appropriate study_id values
train_df <- model_slim %>% filter(study_id %in% train_ids)
valid_df <- model_slim %>% filter(study_id %in% valid_ids)
test_df  <- model_slim %>% filter(study_id %in% test_ids)

# Set seed again
set.seed(2025)

# Ensure levels of sbi outcome are defined appropriately as factors
train_df$SBI <- factor(ifelse(train_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))
valid_df$SBI <- factor(ifelse(valid_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))
test_df$SBI <- factor(ifelse(test_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))

# Make sure character inputs are changed to factors
train_df  <- train_df  %>% mutate(across(where(is.character), as.factor))
valid_df <- valid_df%>% mutate(across(where(is.character), as.factor))
test_df<- test_df%>% mutate(across(where(is.character), as.factor))



# Train model
set.seed(2025)

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 1)  # leave 1 core free

# Work at encounter level to assign folds for 5-fold cross validation within training set
encounter_tbl <- admission_df %>% dplyr::select(study_id, SBI_enc)

# Grouped + stratified 5-fold CV (stratified by sbi outcome)
v <- 5
cv <- rsample::group_vfold_cv(encounter_tbl, v = v, group = study_id, strata = SBI_enc)

# Convert to a simple mapping: study_id -> fold_id (1..v) for the assessment (held-out) side
fold_map <- map2_df(seq_len(v), cv$splits, ~{
  tibble(study_id = assessment(.y)$study_id, fold = .x)
})

# # Save the fold map
# write.csv(fold_map, "/phi/sbi/sbi_blake/cv_fold_map_by_studyid_10_7_25.csv", row.names = FALSE)
# write.csv(fold_map, "/phi/sbi/sbi_blake/cv_fold_map_by_studyid_10_15_25.csv", row.names = FALSE)

# Build caret-style indices from the mapping (for the training set only)
train_fold_map <- fold_map %>%
  filter(study_id %in% train_ids)

# For caret::trainControl(index=...), each element is the row indices used for *training* in that fold.
# We'll include only rows from train_df whose study_id is NOT in that fold's assessment set.
folds <- vector("list", v)
names(folds) <- paste0("Fold", seq_len(v))

for (k in seq_len(v)) {
  heldout_ids <- train_fold_map %>% filter(fold == k) %>% pull(study_id)
  train_rows_k <- which(!train_df$study_id %in% heldout_ids)
  folds[[k]] <- train_rows_k
}

# provide indexOut = held-out rows per fold for transparency
folds_out <- vector("list", v)
names(folds_out) <- paste0("Fold", seq_len(v))
for (k in seq_len(v)) {
  heldout_ids <- train_fold_map %>% filter(fold == k) %>% pull(study_id)
  folds_out[[k]] <- which(train_df$study_id %in% heldout_ids)
}

# Instruct caret to use these folds
ctrl <- trainControl(
  method          = "cv",
  number          = v,              # ignored when index provided
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter     = TRUE,
  index           = folds,
  indexOut        = folds_out,      # optional but still nice to keep
  savePredictions = "final"         # gives out-of-fold preds in rf_tune$pred
)

# Encounter-level split ledger to enable reproducibility
split_ledger <- tibble(
  study_id = all_ids,
  split = case_when(
    study_id %in% train_ids ~ "train",
    study_id %in% valid_ids ~ "valid",
    study_id %in% test_ids  ~ "test",
    TRUE ~ "unknown"
  )
) %>%
  left_join(train_fold_map, by = "study_id") %>%
  mutate(fold = if_else(split == "train", fold, NA_integer_))

# Create split-level summary table with encounter counts, unique patients, and SBI prevalence
encounter_patient_map <- pros_one_model %>%
  dplyr::select(study_id, pat_mrn_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(
    mrn = dplyr::first(pat_mrn_id[!is.na(pat_mrn_id)], default = NA_character_),
    .groups = "drop"
  )

split_summary_table <- split_ledger %>%
  dplyr::select(study_id, split) %>%
  dplyr::left_join(
    admission_df %>% dplyr::select(study_id, sbi_present),
    by = "study_id"
  ) %>%
  dplyr::left_join(encounter_patient_map, by = "study_id") %>%
  dplyr::mutate(
    split = dplyr::recode(
      split,
      train = "training set",
      valid = "validation set",
      test = "test set"
    )
  ) %>%
  dplyr::group_by(split) %>%
  dplyr::summarise(
    n_encounters = dplyr::n_distinct(study_id),
    n_unique_patients = dplyr::n_distinct(mrn, na.rm = TRUE),
    n_encounters_with_sbi = sum(sbi_present == 1, na.rm = TRUE),
    pct_encounters_with_sbi = n_encounters_with_sbi / n_encounters,
    .groups = "drop"
  )

split_summary_table

# # Store csv file with train, test, and validation labels as well as fold labels for training set rows
# write.csv(split_ledger, "/phi/sbi/sbi_blake/split_ledger_studyid_10_7_25.csv", row.names = FALSE)
# write.csv(split_ledger, "/phi/sbi/sbi_blake/split_ledger_studyid_10_15_25.csv", row.names = FALSE)

# Only using the previously identified optimal values for the random forest hyperparameters
rf_grid <- expand.grid(
  mtry            = c(23),
  splitrule       = c("extratrees"),
  min.node.size   = c(3)
)

# Train the model
rf_tune <- train(SBI ~ ., tuneGrid = rf_grid, data = train_df[, c(predictors, "SBI")] , method = "ranger", na.action = na.omit, importance = 'impurity',
                 verbose = TRUE,
                 respect.unordered.factors = TRUE, num.trees = 1000,
                 trControl = ctrl)

# # Extra code to explore hyperparameters
# other_rf <- train(SBI ~ ., data = train_df[, c(predictors, "SBI")] , method = "ranger", na.action = na.omit, importance = 'none',
#                  respect.unordered.factors = TRUE, num.trees = 1000,
#                  trControl = ctrl)
# # Only using the previously identified optimal values for the random forest hyperparameters
# other_grid <- expand.grid(
#   mtry            = c(17:28),
#   splitrule       = c("extratrees", "gini"),
#   min.node.size   = c(1, 3, 5, 10)
# )
#
# other_tune_rf <- train(SBI ~ ., tuneGrid = other_grid, data = train_df[, c(predictors, "SBI")] , method = "ranger", na.action = na.omit, importance = 'impurity',
#                        respect.unordered.factors = TRUE, num.trees = 1000,
#                        trControl = ctrl)




### Optional code to plot top 40 predictors with scaled variable importance
# 1. pull raw importance from the caret model
imp_tbl <- varImp(rf_tune, scale = FALSE)$importance %>%
  tibble::rownames_to_column("predictor") %>%   # keep predictor names
  rename(raw_imp = Overall)

# 2. linearly rescale to 0–100 (100 = most important)
imp_tbl <- imp_tbl %>%
  mutate(scaled_imp = 100 * raw_imp / max(raw_imp, na.rm = TRUE))

# 3. select the top 40 predictors
top40 <- imp_tbl %>%
  arrange(desc(scaled_imp)) %>%
  slice_head(n = 10)

# Change names
top40 %>%
  mutate(predictor = predictor %>%
           str_replace_all("_", " ") %>%    # replace underscores with spaces
           str_to_upper()) %>%              # make all uppercase
  ggplot(aes(x = reorder(predictor, scaled_imp), y = scaled_imp)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Random-Forest Variable Importance (Top 10 Features)",
    x = NULL,
    y = "Scaled Importance (0–100)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),   # larger, bold title
    axis.title.x = element_text(size = 14, face = "bold"), # larger, bold x-axis title
    axis.title.y = element_text(size = 14, face = "bold"), # larger, bold y-axis title
    axis.text = element_text(size = 11)                    # slightly larger tick labels
  )



# Apply to training and test sets
rf_pred_prob_train <- predict(rf_tune, train_df %>% dplyr::select(-study_id, -sbi_present, -abx_exp, -rowid, -SBI), type = "prob")[, "pos"]
rf_pred_prob_val <- predict(rf_tune, valid_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]
rf_pred_prob_test <- predict(rf_tune, test_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]

# Determine AUROC using SBI-negative as the positive/case class.
# Convert sbi_present safely to numeric (0 = no SBI, 1 = SBI present) so subtraction is valid.
to_binary_numeric <- function(x) as.integer(as.character(x))
valid_sbi_neg_case <- 1 - to_binary_numeric(valid_df$sbi_present)
test_sbi_neg_case <- 1 - to_binary_numeric(test_df$sbi_present)

colAUC(1 - rf_pred_prob_val, valid_sbi_neg_case, plotROC = TRUE)
colAUC(1 - rf_pred_prob_test, test_sbi_neg_case, plotROC = TRUE)

# Store these datasets for train, validate, test, along with model predictions
## --- 1) Build data frame for TRAIN set ---

rf_train_df <- train_df %>%
  dplyr::select(study_id, sbi_present, SBI, dplyr::all_of(predictors)) %>%
  dplyr::mutate(rf_prob = rf_pred_prob_train)

## --- 2) Build data frame for VALIDATION set ---

rf_valid_df <- valid_df %>%
  dplyr::select(study_id, sbi_present, SBI, hours_since_picu_adm, dplyr::all_of(predictors)) %>%
  dplyr::mutate(rf_prob = rf_pred_prob_val)

## --- 3) Build data frame for TEST set ---

rf_test_df <- test_df %>%
  dplyr::select(study_id, sbi_present, SBI, hours_since_picu_adm, dplyr::all_of(predictors)) %>%
  dplyr::mutate(rf_prob = rf_pred_prob_test)


# Create calibration plot
# -----------------------------
# Calibration plot on test set
# Predicted probability = Pr(SBI present)
# Observed outcome = SBI present
# -----------------------------

calib_test_df <- rf_test_df %>%
  dplyr::filter(
    !is.na(rf_prob),
    !is.na(sbi_present)
  ) %>%
  dplyr::mutate(
    sbi_present_num = as.integer(as.character(sbi_present))
  )

# Decile-based calibration bins
calib_plot_df <- calib_test_df %>%
  dplyr::mutate(
    pred_bin = dplyr::ntile(rf_prob, 10)
  ) %>%
  dplyr::group_by(pred_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_predicted_risk = mean(rf_prob, na.rm = TRUE),
    observed_sbi_rate = mean(sbi_present_num == 1, na.rm = TRUE),
    .groups = "drop"
  )

p_calibration <- ggplot2::ggplot(
  calib_plot_df,
  ggplot2::aes(x = mean_predicted_risk, y = observed_sbi_rate)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    color = "blue3",
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    color = "blue3",
    size = 3
  ) +
  # ggplot2::geom_text(
  #   ggplot2::aes(label = paste0("n=", n)),
  #   vjust = -0.8,
  #   size = 3.5,
  #   color = "blue3"
  # ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Calibration Plot for New Prospective Random Forest Model",
    subtitle = "Test set; points represent deciles of predicted SBI risk",
    x = "Mean predicted probability of SBI",
    y = "Observed SBI rate"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

p_calibration


######### Now plot NPV, AUROC, and AUPRC by hour #############
# Packages
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(pROC)
library(PRROC)

# -----------------------------
# Settings
# -----------------------------
threshold <- 0.12
n_boot <- 1000
set.seed(123)

# -----------------------------
# Optional preprocessing / duplicate check
# -----------------------------
rf_hour_df <- rf_test_df %>%
  dplyr::filter(hours_since_picu_adm %in% 1:24)

# Overall patient-timepoint pair performance at 0.12 threshold
threshold <- 0.12

overall_sbi_neg_timepoint_summary <- rf_hour_df %>%
  dplyr::filter(
    !is.na(sbi_present),
    !is.na(rf_prob)
  ) %>%
  dplyr::mutate(
    pred_sbi_negative = rf_prob <= threshold,
    actual_sbi_negative = sbi_present == 0,
    correctly_identified_sbi_negative = pred_sbi_negative & actual_sbi_negative
  ) %>%
  dplyr::summarise(
    n_patient_timepoint_pairs = dplyr::n(),
    n_actual_sbi_negative_timepoint_pairs = sum(actual_sbi_negative, na.rm = TRUE),
    n_predicted_sbi_negative_timepoint_pairs = sum(pred_sbi_negative, na.rm = TRUE),
    n_correctly_identified_sbi_negative_timepoint_pairs = sum(correctly_identified_sbi_negative, na.rm = TRUE),
    pct_actual_sbi_negative_timepoints_identified =
      n_correctly_identified_sbi_negative_timepoint_pairs / n_actual_sbi_negative_timepoint_pairs,
    npv_among_predicted_sbi_negative_timepoints =
      n_correctly_identified_sbi_negative_timepoint_pairs / n_predicted_sbi_negative_timepoint_pairs
  )

overall_sbi_neg_timepoint_summary

dup_check <- rf_hour_df %>%
  dplyr::count(study_id, hours_since_picu_adm, name = "n_rows") %>%
  dplyr::filter(n_rows > 1)

if (nrow(dup_check) > 0) {
  warning("Some study_id + hours_since_picu_adm combinations had >1 row. Keeping the first row per patient-hour.")

  rf_hour_df <- rf_hour_df %>%
    dplyr::group_by(study_id, hours_since_picu_adm) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# -----------------------------
# Metric calculation for one dataset
# -----------------------------
calc_metrics_once <- function(dat, threshold = 0.12) {

  dat <- dat %>%
    dplyr::filter(!is.na(sbi_present), !is.na(rf_prob))

  if (nrow(dat) == 0) {
    return(
      tibble::tibble(
        tn = NA_integer_,
        fn = NA_integer_,
        n_negative_pred = NA_integer_,
        n_actual_negative = NA_integer_,
        pct_sbi_neg_identified = NA_real_,
        npv = NA_real_,
        auroc = NA_real_,
        auprc = NA_real_
      )
    )
  }

  y <- dat$sbi_present
  p <- dat$rf_prob
  y_case <- 1 - y
  p_case <- 1 - p

  pred_negative <- p <= threshold

  tn <- sum(pred_negative & y == 0, na.rm = TRUE)
  fn <- sum(pred_negative & y == 1, na.rm = TRUE)
  n_negative_pred <- sum(pred_negative, na.rm = TRUE)

  n_actual_negative <- sum(y == 0, na.rm = TRUE)

  pct_sbi_neg_identified <- if (n_actual_negative > 0) {
    tn / n_actual_negative
  } else {
    NA_real_
  }

  npv <- if (n_negative_pred > 0) {
    tn / n_negative_pred
  } else {
    NA_real_
  }

  auroc <- if (length(unique(y_case)) < 2) {
    NA_real_
  } else {
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = y_case,
          predictor = p_case,
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        )
      )
    )
  }

  auprc <- if (sum(y_case == 1, na.rm = TRUE) == 0 || sum(y_case == 0, na.rm = TRUE) == 0) {
    NA_real_
  } else {
    PRROC::pr.curve(
      scores.class0 = p_case[y_case == 1],
      scores.class1 = p_case[y_case == 0],
      curve = FALSE
    )$auc.integral
  }

  tibble::tibble(
    tn = tn,
    fn = fn,
    n_negative_pred = n_negative_pred,
    n_actual_negative = n_actual_negative,
    pct_sbi_neg_identified = pct_sbi_neg_identified,
    npv = npv,
    auroc = auroc,
    auprc = auprc
  )
}

# -----------------------------
# Bootstrap CI helper
# -----------------------------
get_boot_ci <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) {
    return(c(NA_real_, NA_real_))
  }

  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

# Split-level performance at threshold = 0.12 with bootstrap CIs
format_metric_ci <- function(est, lower, upper, digits = 3) {
  if (is.na(est) || is.na(lower) || is.na(upper)) {
    return(NA_character_)
  }
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), est, lower, upper)
}

calc_split_metrics_boot <- function(dat, threshold = 0.12, n_boot = 1000) {
  dat <- dat %>%
    dplyr::filter(!is.na(sbi_present), !is.na(rf_prob))

  if (nrow(dat) == 0) {
    return(tibble::tibble(
      auroc = NA_real_, auroc_lower = NA_real_, auroc_upper = NA_real_,
      auprc = NA_real_, auprc_lower = NA_real_, auprc_upper = NA_real_,
      npv = NA_real_, npv_lower = NA_real_, npv_upper = NA_real_
    ))
  }

  point_est <- calc_metrics_once(dat, threshold = threshold)

  boot_res <- replicate(
    n = n_boot,
    expr = {
      boot_idx <- sample.int(n = nrow(dat), size = nrow(dat), replace = TRUE)
      boot_dat <- dat[boot_idx, , drop = FALSE]
      boot_met <- calc_metrics_once(boot_dat, threshold = threshold)
      c(auroc = boot_met$auroc, auprc = boot_met$auprc, npv = boot_met$npv)
    },
    simplify = "matrix"
  )

  auroc_ci <- get_boot_ci(boot_res["auroc", ])
  auprc_ci <- get_boot_ci(boot_res["auprc", ])
  npv_ci <- get_boot_ci(boot_res["npv", ])

  tibble::tibble(
    auroc = point_est$auroc,
    auroc_lower = auroc_ci[1],
    auroc_upper = auroc_ci[2],
    auprc = point_est$auprc,
    auprc_lower = auprc_ci[1],
    auprc_upper = auprc_ci[2],
    npv = point_est$npv,
    npv_lower = npv_ci[1],
    npv_upper = npv_ci[2]
  )
}

split_performance_table <- dplyr::bind_rows(
  "training set" = calc_split_metrics_boot(rf_train_df, threshold = 0.12, n_boot = 1000),
  "validation set" = calc_split_metrics_boot(rf_valid_df, threshold = 0.12, n_boot = 1000),
  "test set" = calc_split_metrics_boot(rf_test_df, threshold = 0.12, n_boot = 1000),
  .id = "split"
) %>%
  dplyr::mutate(
    `AUROC (95% CI)` = format_metric_ci(auroc, auroc_lower, auroc_upper),
    `AUPRC (95% CI)` = format_metric_ci(auprc, auprc_lower, auprc_upper),
    `NPV (95% CI)` = format_metric_ci(npv, npv_lower, npv_upper)
  ) %>%
  dplyr::select(split, `AUROC (95% CI)`, `AUPRC (95% CI)`, `NPV (95% CI)`)

split_summary_table <- split_summary_table %>%
  dplyr::left_join(split_performance_table, by = "split")

split_summary_table



# -----------------------------
# Compute point estimates + bootstrap CIs for one hour
# -----------------------------
calc_hour_metrics_boot <- function(dat, threshold = 0.12, n_boot = 1000) {

  dat <- dat %>%
    dplyr::filter(!is.na(sbi_present), !is.na(rf_prob))

  if (nrow(dat) == 0) {
    return(
      tibble::tibble(
        n_obs = 0,
        n_patients = 0,
        tn = NA_integer_,
        fn = NA_integer_,
        n_negative_pred = NA_integer_,
        n_actual_negative = NA_integer_,
        pct_sbi_neg_identified = NA_real_,
        npv = NA_real_,
        npv_lower = NA_real_,
        npv_upper = NA_real_,
        auroc = NA_real_,
        auroc_lower = NA_real_,
        auroc_upper = NA_real_,
        auprc = NA_real_,
        auprc_lower = NA_real_,
        auprc_upper = NA_real_
      )
    )
  }

  point_est <- calc_metrics_once(dat, threshold = threshold)

  boot_res <- replicate(
    n = n_boot,
    expr = {
      boot_idx <- sample.int(n = nrow(dat), size = nrow(dat), replace = TRUE)
      boot_dat <- dat[boot_idx, , drop = FALSE]
      boot_met <- calc_metrics_once(boot_dat, threshold = threshold)
      c(
        npv = boot_met$npv,
        auroc = boot_met$auroc,
        auprc = boot_met$auprc
      )
    },
    simplify = "matrix"
  )

  npv_ci <- get_boot_ci(boot_res["npv", ])
  auroc_ci <- get_boot_ci(boot_res["auroc", ])
  auprc_ci <- get_boot_ci(boot_res["auprc", ])

  tibble::tibble(
    n_obs = nrow(dat),
    n_patients = dplyr::n_distinct(dat$study_id),
    tn = point_est$tn,
    fn = point_est$fn,
    n_negative_pred = point_est$n_negative_pred,
    n_actual_negative = point_est$n_actual_negative,
    pct_sbi_neg_identified = point_est$pct_sbi_neg_identified,
    npv = point_est$npv,
    npv_lower = npv_ci[1],
    npv_upper = npv_ci[2],
    auroc = point_est$auroc,
    auroc_lower = auroc_ci[1],
    auroc_upper = auroc_ci[2],
    auprc = point_est$auprc,
    auprc_lower = auprc_ci[1],
    auprc_upper = auprc_ci[2]
  )
}

# -----------------------------
# Run for each hour 1 through 24
# -----------------------------
rf_hourly_metrics_boot <- purrr::map_dfr(
  1:24,
  function(hr) {
    dat_hr <- rf_hour_df %>%
      dplyr::filter(hours_since_picu_adm == hr)

    calc_hour_metrics_boot(
      dat = dat_hr,
      threshold = threshold,
      n_boot = n_boot
    ) %>%
      dplyr::mutate(hours_since_picu_adm = hr, .before = 1)
  }
) %>%
  dplyr::mutate(
    auroc = round(auroc, 2),
    auroc_lower = round(auroc_lower, 2),
    auroc_upper = round(auroc_upper, 2),
    auprc = round(auprc, 2),
    auprc_lower = round(auprc_lower, 2),
    auprc_upper = round(auprc_upper, 2),
    npv_label = dplyr::if_else(
      !is.na(n_negative_pred) & n_negative_pred > 0,
      paste0(tn, "/", n_negative_pred),
      NA_character_
    ),
    pct_sbi_neg_identified_label = dplyr::if_else(
      !is.na(pct_sbi_neg_identified),
      paste0(round(100 * pct_sbi_neg_identified), "%"),
      NA_character_
    ),
    auroc_label = dplyr::if_else(
      !is.na(auroc),
      sprintf("%.2f", auroc),
      NA_character_
    ),
    auprc_label = dplyr::if_else(
      !is.na(auprc),
      sprintf("%.2f", auprc),
      NA_character_
    )
  )

# View results table
rf_hourly_metrics_boot

# Compute SBI-negative prevalence for AUPRC plot
# IMPORTANT:
# AUPRC is being calculated with SBI-negative as the case/event class.
# Therefore, the no-skill / prevalence reference line should be:
#   number of SBI-negative encounters / total encounters
# not SBI-positive prevalence.

# Helper to safely coerce sbi_present to numeric 0/1
# Expected coding:
#   0 = SBI negative
#   1 = SBI positive
to_sbi_numeric <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x_clean <- stringr::str_to_lower(stringr::str_trim(x))

    return(dplyr::case_when(
      x_clean %in% c("0", "no", "neg", "negative", "sbi-", "sbi_negative", "false") ~ 0L,
      x_clean %in% c("1", "yes", "pos", "positive", "sbi+", "sbi_positive", "true") ~ 1L,
      TRUE ~ NA_integer_
    ))
  }

  as.integer(x)
}

# Check that sbi_present is constant within each patient/encounter
sbi_consistency_check <- rf_hour_df %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present)
  ) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(
    n_outcome_values = dplyr::n_distinct(sbi_present_num[!is.na(sbi_present_num)]),
    .groups = "drop"
  )

if (any(sbi_consistency_check$n_outcome_values > 1)) {
  warning("Some study_id values have more than one sbi_present value across hours.")
}

# One row per patient/encounter for prevalence calculation
patient_level_outcome <- rf_hour_df %>%
  dplyr::mutate(
    sbi_present_num = to_sbi_numeric(sbi_present)
  ) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(
    sbi_present_num = dplyr::first(sbi_present_num[!is.na(sbi_present_num)]),
    .groups = "drop"
  )

# Explicit SBI-negative prevalence
sbi_neg_prevalence_df <- patient_level_outcome %>%
  dplyr::summarise(
    n_encounters = dplyr::n(),
    n_sbi_negative = sum(sbi_present_num == 0, na.rm = TRUE),
    n_sbi_positive = sum(sbi_present_num == 1, na.rm = TRUE),
    sbi_neg_prevalence = n_sbi_negative / n_encounters,
    sbi_pos_prevalence = n_sbi_positive / n_encounters
  )

sbi_neg_prevalence <- sbi_neg_prevalence_df$sbi_neg_prevalence

sbi_neg_prevalence_label <- paste0(
  "SBI-negative prevalence = ",
  scales::percent(sbi_neg_prevalence, accuracy = 1)
)

# Optional QC printout
sbi_neg_prevalence_df



## Plot NPV ##
p_npv <- ggplot2::ggplot(
  rf_hourly_metrics_boot,
  ggplot2::aes(x = hours_since_picu_adm)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = npv_lower, ymax = npv_upper),
    fill = "lightblue",
    alpha = 0.4,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = npv),
    color = "blue3",
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = npv),
    color = "blue3",
    size = 2.5,
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = npv, label = npv_label),
    color = "blue3",
    vjust = -0.8,
    size = 3.5,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = pct_sbi_neg_identified),
    color = "firebrick3",
    linewidth = 0.9,
    linetype = "dashed",
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = pct_sbi_neg_identified),
    color = "firebrick3",
    size = 2.3,
    shape = 17,
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = pct_sbi_neg_identified,
      label = pct_sbi_neg_identified_label
    ),
    color = "firebrick3",
    vjust = 1.5,
    size = 3.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_x_continuous(breaks = 1:24) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.1),
    expand = ggplot2::expansion(mult = c(0.01, 0.10)),
    name = "NPV",
    sec.axis = ggplot2::sec_axis(
      transform = ~ . * 100,
      name = "% of SBI-negative patients correctly identified"
    )
  ) +
  ggplot2::labs(
    title = "Negative Predictive Value by PICU Hour",
    x = "Hours Since PICU Admission"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.title.x = ggplot2::element_text(face = "bold"),
    axis.title.y.left = ggplot2::element_text(face = "bold", color = "blue3"),
    axis.title.y.right = ggplot2::element_text(face = "bold", color = "firebrick3"),
    axis.text.y.left = ggplot2::element_text(color = "blue3"),
    axis.text.y.right = ggplot2::element_text(color = "firebrick3"),
    panel.grid.minor = ggplot2::element_blank()
  )

p_npv

## Plot AUROC ##
p_auroc <- ggplot2::ggplot(
  rf_hourly_metrics_boot,
  ggplot2::aes(x = hours_since_picu_adm, y = auroc)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = auroc_lower, ymax = auroc_upper),
    fill = "lightblue",
    alpha = 0.4,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    color = "blue3",
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    color = "blue3",
    size = 2.5,
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = auroc_label),
    color = "blue3",
    vjust = -0.8,
    size = 3.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_x_continuous(breaks = 1:24) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.1),
    expand = ggplot2::expansion(mult = c(0.01, 0.08))
  ) +
  ggplot2::labs(
    title = "AUROC by PICU Hour",
    x = "Hours Since PICU Admission",
    y = "AUROC"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

p_auroc


## AUPRC plot ##
p_auprc <- ggplot2::ggplot(
  rf_hourly_metrics_boot,
  ggplot2::aes(x = hours_since_picu_adm, y = auprc)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = auprc_lower, ymax = auprc_upper),
    fill = "lightblue",
    alpha = 0.4,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    color = "blue3",
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    color = "blue3",
    size = 2.5,
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = auprc_label),
    color = "blue3",
    vjust = -0.8,
    size = 3.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_x_continuous(breaks = 1:24) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.1),
    expand = ggplot2::expansion(mult = c(0.01, 0.08))
  ) +
  ggplot2::labs(
    title = "AUPRC by PICU Hour",
    x = "Hours Since PICU Admission",
    y = "AUPRC"
  ) +

  ggplot2::geom_hline(
    yintercept = sbi_neg_prevalence,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  ggplot2::annotate(
    "text",
    x = 24,
    y = sbi_neg_prevalence,
    label = sbi_neg_prevalence_label,
    hjust = 1,
    vjust = -0.6,
    size = 3.8,
    color = "black",
    fontface = "bold"
  ) +

  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

p_auprc
