# Script to evaluate new model that doesn't leverage unit specific origin locations as input variables #

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

# At this point there are 69 predictors, will remove the preicu variables for this sensitivity analysis
remove_pre_icu <- c("preicu_9thfloor", "preicu_er", "preicu_nocer", "preicu_nonchco", "preicu_or")

predictors <- predictors[!predictors %in% remove_pre_icu]

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
# write.csv(fold_map, file.path(sbi_blake_phi_path, "cv_fold_map_by_studyid_10_7_25.csv"), row.names = FALSE)
# write.csv(fold_map, file.path(sbi_blake_phi_path, "cv_fold_map_by_studyid_10_15_25.csv"), row.names = FALSE)

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
encounter_patient_map <- pros_one_model %>% mutate(pat_mrn_id = as.character(pat_mrn_id)) %>%
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

split_summary_table <- split_summary_table %>%
  dplyr::mutate(
    split = factor(
      split,
      levels = c("training set", "validation set", "test set")
    ),
    pct_encounters_with_sbi = round(100 * pct_encounters_with_sbi, 1)
  ) %>%
  dplyr::arrange(split) %>%
  dplyr::mutate(
    split = as.character(split)
  )



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
# write.csv(fold_map, file.path(sbi_blake_phi_path, "cv_fold_map_by_studyid_10_7_25.csv"), row.names = FALSE)
# write.csv(fold_map, file.path(sbi_blake_phi_path, "cv_fold_map_by_studyid_10_15_25.csv"), row.names = FALSE)

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
encounter_patient_map <- pros_one_model %>% mutate(pat_mrn_id = as.character(pat_mrn_id)) %>%
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

split_summary_table <- split_summary_table %>%
  dplyr::mutate(
    split = factor(
      split,
      levels = c("training set", "validation set", "test set")
    ),
    pct_encounters_with_sbi = round(100 * pct_encounters_with_sbi, 1)
  ) %>%
  dplyr::arrange(split) %>%
  dplyr::mutate(
    split = as.character(split)
  )

# Explore hyperparameters
other_rf <- train(SBI ~ ., data = train_df[, c(predictors, "SBI")] , method = "ranger", na.action = na.omit, importance = 'none',
                 respect.unordered.factors = TRUE, num.trees = 1000,
                 trControl = ctrl)
