#### Script to try and design different models using retro data ####

### Script to create the retrospective models using the correct dataset and then apply to the prospective inputs for
# updated model probabilities


model_data_df <- read.fst(path = "~/sbi_blake/jan_25_23_model_data_df.fst")
model_glm_orig <- model_data_df

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

rf_df <- model_data_df %>% left_join(case_id_and_picu_adm_date_time, by = "case_id")  %>% filter(picu_adm_date_time < t_beaker) %>% filter(abx_exp == "0") %>%
  dplyr::select(-all_of(c("picu_adm_date_time"))) #n = 8,657

bad_cols <- c("sbi_pneumonia", "sbi", "blood", "abd", "cns", "genital", "nos", "resp", "stool", "tissue", "urine", "sbi_cx_neg_sepsis", "virus_24", "abx_exp", "death_date")
rf_df <- rf_df %>% dplyr::select(-all_of(bad_cols))

# Optional code to only include the top 40 variables identified later
au_40_list <- c(
  "mean_sat",	"mean_fio2",	"min_fio2",	"mean_temp",	"slope_sat",	"slope_rr",	"min_sbp",	"slope_hr",	"median_temp",	"median_fio2",
  "mean_sbp",	"median_sbp",	"last_sbp",	"max_temp",	"age_years",	"last_fio2",	"min_temp",	"max_sbp",	"slope_dbp",	"min_dbp",	"min_hr",
  "max_dbp",	"slope_sbp",	"last_hr",	"last_dbp",	"los_before_icu_days",	"mean_dbp",	"number_fio2",	"median_hr",	"mean_hr",	"max_fio2",
  "median_sat",	"median_dbp",	"max_hr",	"last_temp",	"slope_temp",	"min_sat",	"slope_fio2",	"min_rr",	"pco2_gas_blood"
)

rf_df <- rf_df %>% dplyr::select(all_of(c("case_id", "sbi_present", au_40_list)))

# #################################################################################################################################################

#Fix classes
rf_df$sbi_present <- factor(
  rf_df$sbi_present,
  levels = c("yes", "no")
)

#################################################################################################################################################

rf_df <- rf_df %>% filter(max_sbp <= 100)

rf_new_names <- rf_df %>%
  dplyr::rename(
    o2sat_mean       = mean_sat,
    fio2_mean        = mean_fio2,
    fio2_min         = min_fio2,
    temp_mean        = mean_temp,
    o2sat_slope      = slope_sat,
    rr_slope         = slope_rr,
    sbp_min          = min_sbp,
    hr_slope         = slope_hr,
    temp_median      = median_temp,
    fio2_median      = median_fio2,
    sbp_mean         = mean_sbp,
    sbp_median       = median_sbp,
    sbp_last         = last_sbp,
    temp_max         = max_temp,
    age              = age_years,
    fio2_last        = last_fio2,
    temp_min         = min_temp,
    sbp_max          = max_sbp,
    dbp_slope        = slope_dbp,
    dbp_min          = min_dbp,
    hr_min           = min_hr,
    dbp_max          = max_dbp,
    sbp_slope        = slope_sbp,
    hr_last          = last_hr,
    dbp_last         = last_dbp,
    los_before_icu_days = los_before_icu_days,
    dbp_mean         = mean_dbp,
    fio2_count       = number_fio2,
    hr_median        = median_hr,
    hr_mean          = mean_hr,
    fio2_max         = max_fio2,
    o2sat_median     = median_sat,
    dbp_median       = median_dbp,
    hr_max           = max_hr,
    temp_last        = last_temp,
    temp_slope       = slope_temp,
    o2sat_min        = min_sat,
    fio2_slope       = slope_fio2,
    rr_min           = min_rr,
    pco2_mean        = pco2_gas_blood
  )


# First divide the rf_df into the train and test data frames by random 80/20 split
set.seed(seed = 100) #32313
train_size <- floor(nrow(rf_new_names) * 0.8)
sample <- sample.int(n = nrow(rf_new_names), size = train_size, replace = FALSE)
train_df <- rf_new_names[sample, ] #6920
test_df <- rf_new_names[-sample, ] # 1731 patients

predictors <- setdiff(names(train_df), c("case_id", "sbi_present"))

form <- as.formula(
  paste("sbi_present ~", paste(predictors, collapse = " + "))
)
#
rf_model_trial <- train(
  form,
  tuneLength = 10,
  data = train_df,
  method = "ranger",
  na.action = na.omit,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  trControl = trainControl(
    method = "cv", number = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE
  )
)

new_mtry <- 5; new_split <- "gini";
#
myGrid <- data.frame(.mtry = rep(c((new_mtry - 4):(new_mtry + 5)), times = 4), .splitrule = rep(c(as.character(new_split)), each = 40), .min.node.size = rep(c(1:4), each = 10))

rf_model <- train(form, tuneGrid = myGrid , data = (train_df), method = "ranger", na.action = na.omit, importance = "impurity",
                                     respect.unordered.factors = TRUE,
                                     trControl = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE, verboseIter = TRUE))

# 3, gini, 1 best parameters with OOF AUROC of 0.715

# Rerun with just these parameters and save out of fold probabilities
myGrid <- data.frame(.mtry = 3, .splitrule = "gini", .min.node.size = 1)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = TRUE,
  savePredictions = "final"   # <-- THIS is the key
)

rf_model <- train(
  form,
  tuneGrid = myGrid,
  data = train_df,
  method = "ranger",
  na.action = na.omit,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  trControl = ctrl
)

# Store training data OOF AUROC which is 0.713
retro_train_auroc_no_abx <- rf_model$results %>%
  dplyr::filter(
    mtry == rf_model$bestTune$mtry,
    min.node.size == rf_model$bestTune$min.node.size
  ) %>%
  dplyr::pull(ROC)


# Determine NPV in training set
# ---- 1) Determine which probability column corresponds to the "event"/positive class ----
# caret stores probs as columns named by class labels. We want P(Y=positive).
get_positive_prob_col <- function(model) {
  lev <- model$levels
  # by default, caret uses the FIRST level as the "event" for twoClassSummary
  # (this is the one used for ROC). So probability col is lev[1].
  lev[1]
}

# ---- 2) Compute NPV (and a few other useful quantities) at a threshold ----
npv_at_threshold <- function(model, threshold = 0.05) {
  pred <- model$pred
  stopifnot(!is.null(pred), nrow(pred) > 0)

  pos_lab <- model$levels[1]
  neg_lab <- model$levels[2]
  prob_col <- get_positive_prob_col(model)

  # If you did repeated CV / multiple tunes, keep bestTune only
  if (!is.null(model$bestTune)) {
    pred <- pred %>% inner_join(model$bestTune, by = names(model$bestTune))
  }

  d <- pred %>%
    transmute(
      truth = obs,
      p_pos = .data[[prob_col]],
      pred_neg = (p_pos <= threshold)  # "rule-out" / predicted negative
    ) %>%
    filter(!is.na(truth), !is.na(p_pos))

  # Confusion elements for "pred_neg" group
  tn <- sum(d$pred_neg & d$truth == neg_lab)
  fn <- sum(d$pred_neg & d$truth == pos_lab)
  n_pred_neg <- sum(d$pred_neg)

  npv <- if (n_pred_neg > 0) tn / (tn + fn) else NA_real_

  tibble(
    threshold = threshold,
    n = nrow(d),
    n_pred_neg = n_pred_neg,
    tn = tn,
    fn = fn,
    npv = npv,
    fn_per_1000 = if (is.finite(npv)) 1000 * (fn / n_pred_neg) else NA_real_
  )
}

# ---- 3) Find largest threshold with OOF NPV meeting target ----
find_largest_threshold_for_npv <- function(model, target_npv = 0.95, threshold_grid = seq(0, 1, by = 0.001)) {
  npv_curve <- purrr::map_dfr(threshold_grid, function(t) npv_at_threshold(model, threshold = t))

  eligible <- npv_curve %>%
    dplyr::filter(!is.na(npv), n_pred_neg > 0, npv >= target_npv) %>%
    dplyr::arrange(dplyr::desc(threshold))

  if (nrow(eligible) == 0) {
    best_available <- npv_curve %>%
      dplyr::filter(!is.na(npv), n_pred_neg > 0) %>%
      dplyr::arrange(dplyr::desc(npv), dplyr::desc(n_pred_neg), threshold) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(selection_rule = "max_npv_fallback")

    return(best_available)
  }

  eligible %>%
    dplyr::slice(1) %>%
    dplyr::mutate(selection_rule = "largest_threshold_meeting_target_npv")
}

# Shows out of fold NPV at a given threshold
npv_at_threshold(rf_model, threshold = 0.050) # 0.958 = NPV


# Now will apply this model to the test set and compute the AUROC and associated confusion matrix
rf_pred_prob <- predict(rf_model, test_df, type = "prob")
retro_test_auroc_no_abx <- as.numeric(
  pROC::auc(
    pROC::roc(
      response  = test_df$sbi_present,
      predictor = rf_pred_prob[["yes"]],
      levels    = c("no", "yes"),
      direction = "<",
      quiet     = TRUE
    )
  )
)

# Determine the largest OOF threshold that still provides NPV >= 0.95
oof_threshold_no_abx <- find_largest_threshold_for_npv(
  rf_model,
  target_npv = 0.95
)
threshold_no_abx <- oof_threshold_no_abx$threshold
retro_train_npv_no_abx <- oof_threshold_no_abx$npv

# Apply the selected OOF-derived threshold to the holdout test set
pred_no_abx <- ifelse(rf_pred_prob[["yes"]] <= threshold_no_abx, "no", "yes")

rf_cm_no_abx_test <- caret::confusionMatrix(
  data = factor(pred_no_abx, levels = c("no", "yes")),
  reference = factor(as.character(test_df$sbi_present), levels = c("no", "yes")),
  mode = "everything",
  positive = "yes"
)

npv_retro_no_abx_test <- unname(rf_cm_no_abx_test$byClass["Neg Pred Value"])

# Store selected test performance summary
retro_test_threshold_selection_no_abx <- tibble::tibble(
  threshold_source = as.character(oof_threshold_no_abx$selection_rule),
  selected_threshold = threshold_no_abx,
  selected_npv_test = npv_retro_no_abx_test
)


###### Exrta code to look at NPV by threshold #######
## NPV by threshold
pred_prob_no_abx <- rf_pred_prob[["yes"]]
truth_no_abx <- as.character(test_df$sbi_present)

thresholds <- seq(0.01, 0.15, by = 0.001)
npv_by_threshold <- purrr::map_dfr(thresholds, function(ocp_rf) {

  pred_class <- ifelse(pred_prob_no_abx >= ocp_rf, "yes", "no")

  tn <- sum(pred_class == "no" & truth_no_abx == "no", na.rm = TRUE)
  fn <- sum(pred_class == "no" & truth_no_abx == "yes", na.rm = TRUE)
  n_pred_negative <- tn + fn

  tibble::tibble(
    threshold = ocp_rf,
    npv = ifelse(n_pred_negative > 0, tn / n_pred_negative, NA_real_),
    n_pred_negative = n_pred_negative,
    n_true_negative = tn,
    n_false_negative = fn
  )
})

View(npv_by_threshold)


# Store AUPRC in test set also
# make a data frame with truth + positive-class prob
test_pr_df <- tibble(
  truth = test_df$sbi_present,
  p_yes = rf_pred_prob[, "yes"]   # <-- positive class prob (adjust if your level name differs)
) %>%
  mutate(truth = factor(truth, levels = c("yes","no")))

retro_test_auprc_no_abx <- yardstick::pr_auc(test_pr_df, truth = truth, p_yes) %>% # 0.51 in test set, prevalence 0.256
  pull(.estimate)


# Store the training set AURPC
# oof predictions saved by caret (requires savePredictions="final")
oof <- rf_model$pred

# AUPRC using PRROC (class0 = positive class scores)
oof <- rf_model$pred %>% mutate(obs = factor(obs, levels = c("yes","no")))
retro_train_auprc_no_abx <- yardstick::pr_auc(
  oof,
  truth = obs,
  yes        # <-- bare column name, no `estimate =`
) %>%
  pull(.estimate) #0.486 AUPRC, prevalence in training set of SBI is 0.269


# Store OOF probs
oof_preds <- rf_model$pred %>%
  dplyr::select(rowIndex, obs, yes) %>%
  dplyr::rename(
    sbi_present = obs,
    oof_prob = yes
  )

train_df_with_oof <- train_df %>%
  dplyr::mutate(rowIndex = dplyr::row_number()) %>%
  dplyr::left_join(oof_preds, by = "rowIndex")

# Fix and store train_df_with_oof
train_df_with_oof <- train_df_with_oof %>% dplyr::select(-sbi_present.y) %>% rename(sbi_present = sbi_present.x)
train_df_with_oof <- train_df_with_oof %>% rename(model_score = oof_prob) %>% dplyr::select(-rowIndex)


# -----------------------------
# Create retro TEST set dataframe with inputs + predicted probability
# -----------------------------

# 1) Get predicted probabilities on the holdout test set (you already do this)
rf_pred_prob <- predict(rf_model, test_df, type = "prob")

# 2) Build a tidy predictions df keyed by row order
test_preds <- tibble(
  rowIndex   = seq_len(nrow(test_df)),                  # 1:nrow(test_df)
  sbi_present = test_df$sbi_present,                    # keep factor yes/no
  model_score = rf_pred_prob[, "yes"]                   # P(SBI = yes)
)

roc_obj <- roc(
  response  = test_preds$sbi_present,
  predictor = test_preds$model_score,
  levels = c("no", "yes"),   # sets "yes" as the positive class
  direction = "<"
)

auc(roc_obj) # AUROC 0.745 in test set

# Now calculate NPV at OOF-selected threshold

test_preds_npv <- test_preds %>%
  mutate(
    pred_class = if_else(model_score <= threshold_no_abx, "no", "yes"),
    pred_class = factor(pred_class, levels = c("no", "yes")),
    sbi_present = factor(sbi_present, levels = c("no", "yes"))
  )

# NPV = TN / (TN + FN)
npv_retro_no_abx <- with(test_preds_npv,
                         sum(pred_class == "no" & sbi_present == "no") /
                           sum(pred_class == "no"))

npv_retro_no_abx # NPV 0.93


# 3) Attach to the test_df inputs
test_df_with_pred <- test_df %>%
  mutate(rowIndex = dplyr::row_number()) %>%
  left_join(test_preds %>% select(rowIndex, model_score), by = "rowIndex") %>%
  select(-rowIndex) %>%
  relocate(case_id, sbi_present, model_score)

# Sanity checks
stopifnot(nrow(test_df_with_pred) == nrow(test_df))
stopifnot(!anyNA(test_df_with_pred$model_score))

glimpse(test_df_with_pred)

rf_df <- test_df_with_pred

# write_csv(x = train_df_with_oof, file = "/home/martinbl/sbi_blake/validation_files/rf_au_40_training_updated_2_4_26.csv")

############################################ Identical code to create the abx_exposed version of the model #################################
rf_df_abx <- model_data_df %>% left_join(case_id_and_picu_adm_date_time, by = "case_id")  %>%
  filter(picu_adm_date_time < t_beaker) %>% filter(abx_exp == "1") %>%
  dplyr::select(-all_of(c("picu_adm_date_time"))) #n = 8,657

bad_cols <- c("sbi_pneumonia", "sbi", "blood", "abd", "cns", "genital", "nos", "resp", "stool", "tissue", "urine",
              "sbi_cx_neg_sepsis", "virus_24", "abx_exp", "death_date")

rf_df_abx <- rf_df_abx %>% dplyr::select(-all_of(bad_cols))

#Fix classes
rf_df_abx$sbi_present <- factor(
  rf_df_abx$sbi_present,
  levels = c("yes", "no")
)

# Include the top 40 variables identified later
ae_40_list <- c(
  "age_years",
  "hematocrit_blood",
  "lactate_blood",
  "last_dbp",
  "last_fio2",
  "last_hr",
  "last_sbp",
  "last_temp",
  "los_before_icu_days",
  "max_dbp",
  "max_hr",
  "max_sbp",
  "max_temp",
  "mean_dbp",
  "mean_fio2",
  "mean_hr",
  "mean_rr",
  "mean_sat",
  "mean_sbp",
  "mean_temp",
  "median_dbp",
  "median_fio2",
  "median_hr",
  "median_rr",
  "median_sbp",
  "median_temp",
  "min_dbp",
  "min_hr",
  "min_rr",
  "min_sat",
  "min_sbp",
  "number_sat",
  "pre_icu",
  "slope_dbp",
  "slope_hr",
  "slope_rr",
  "slope_sat",
  "slope_sbp",
  "slope_temp"
)

rf_df_abx <- rf_df_abx %>% dplyr::select(all_of(c("case_id", "sbi_present", ae_40_list)))

#Fix classes
rf_df_abx$sbi_present <- as.factor(rf_df_abx$sbi_present)


rf_df_abx <- rf_df_abx %>% filter(max_sbp <= 100)

# One hot encode the pre_icu column
rf_df_abx$pre_icu <- factor(rf_df_abx$pre_icu)

# One-hot encode
pre_icu_ohe <- model.matrix(~ pre_icu - 1, data = rf_df_abx)

# Clean column names
colnames(pre_icu_ohe) <- gsub("^pre_icu", "pre_icu_", colnames(pre_icu_ohe))

# Bind back to dataframe (optionally drop original column)
rf_df_abx <- rf_df_abx %>%
  dplyr::select(-pre_icu) %>%
  dplyr::bind_cols(as.data.frame(pre_icu_ohe))


rf_new_names_abx <- rf_df_abx %>%
  dplyr::rename(
    o2sat_mean       = mean_sat,
    lactate_max = lactate_blood,
    rr_mean = mean_rr,
    preicu_6thfloor = pre_icu_6th_floor,
    preicu_8thfloor = pre_icu_8th_floor,
    preicu_9thfloor = pre_icu_9th_floor,
    preicu_er = pre_icu_er,
    preicu_hem_onc_bmt = pre_icu_hem_onc_bmt,
    preicu_nocer = pre_icu_noc_er,
    preicu_nocinpatient = pre_icu_noc_inpatient,
    preicu_nonchco = `pre_icu_Non-CHCO Facility`,
    preicu_or = pre_icu_operating_room,
    preicu_othericu = pre_icu_other_CHCO_ICU,
    preicu_outpatient = pre_icu_outpatient,
    preicu_procedure_center = pre_icu_procedure_center,
    sbp_min = min_sbp,
    rr_median = median_rr,
    los_before_icu_days= los_before_icu_days,
    fio2_mean= mean_fio2,
    rr_slope = slope_rr,
    o2sat_min = min_sat,
    sbp_mean = mean_sbp,
    hr_slope = slope_hr,
    sbp_median = median_sbp,
    temp_max = max_temp,
    temp_slope = slope_temp,
    o2sat_slope = slope_sat,
    fio2_median = median_fio2,
    dbp_min = min_dbp,
    sbp_slope = slope_sbp,
    temp_mean = mean_temp,
    dbp_slope = slope_dbp,
    sbp_last = last_sbp,
    age = age_years,
    hr_min = min_hr,
    o2sat_count= number_sat,
    hr_max = max_hr,
    hr_mean = mean_hr,
    temp_median = median_temp,
    dbp_median = median_dbp,
    hr_median= median_hr,
    dbp_mean = mean_dbp,
    hr_last = last_hr,
    rr_min = min_rr,
    sbp_max = max_sbp,
    dbp_last = last_dbp,
    temp_last = last_temp
  )

# Divide the rf_df_abx into the train and test data frames by random 80/20 split
set.seed(seed = 100) #32313
train_size <- floor(nrow(rf_new_names_abx) * 0.8)
sample <- sample.int(n = nrow(rf_new_names_abx), size = train_size, replace = FALSE)
train_df_abx <- rf_new_names_abx[sample, ] #3302
test_df_abx <- rf_new_names_abx[-sample, ] # 826 patients

predictors_abx <- setdiff(names(train_df_abx), c("case_id", "sbi_present"))

formabx <- as.formula(
  paste("sbi_present ~", paste(predictors_abx, collapse = " + "))
)

rf_model_trial <- train(
  formabx,
  tuneLength = 10,
  data = train_df_abx,
  method = "ranger",
  na.action = na.omit,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  trControl = trainControl(
    method = "cv", number = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE
  )
)

new_mtry <- 18; new_split <- "gini";

myGrid <- data.frame(.mtry = rep(c((new_mtry - 4):(new_mtry + 5)), times = 4), .splitrule = rep(c(as.character(new_split)), each = 40), .min.node.size = rep(c(1:4), each = 10))
# 23, gini, 3 best parameters with OOF AUROC of 0.76

# Rerun with just these parameters and save out of fold probabilities
myGrid <- data.frame(.mtry = 23, .splitrule = "gini", .min.node.size = 3)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = TRUE,
  savePredictions = "final"   # <-- THIS is the key
)

rf_model_abx <- train(
  formabx,
  tuneGrid = myGrid,
  data = train_df_abx,
  method = "ranger",
  na.action = na.omit,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  trControl = ctrl
)

# Store training data OOF AUROC which is 0.713
retro_train_auroc_yes_abx <- rf_model_abx$results %>%
  dplyr::filter(
    mtry == rf_model_abx$bestTune$mtry,
    min.node.size == rf_model_abx$bestTune$min.node.size
  ) %>%
  dplyr::pull(ROC)

oof <- rf_model_abx$pred

# Keep only the row for the bestTune (important if you ever used tuneLength / multiple grid rows)
oof <- oof %>%
  semi_join(rf_model_abx$bestTune, by = names(rf_model_abx$bestTune))

# Identify largest OOF threshold that still gives NPV >= 0.95
oof_threshold_yes_abx <- find_largest_threshold_for_npv(
  rf_model_abx,
  target_npv = 0.95
)
threshold_yes_abx <- oof_threshold_yes_abx$threshold
retro_train_npv_yes_abx <- oof_threshold_yes_abx$npv


# Now will apply this model to the training set and compute the AUROC and associated confusion matrix
rf_pred_prob_abx <- predict(rf_model_abx, test_df_abx, type = "prob")

### Look at NPV by threshold in test set
# Predicted probability of SBI = yes for the antibiotic-exposed cohort
pred_prob_abx <- rf_pred_prob_abx[["yes"]]

# True outcome as character, preserving existing yes/no labels
truth_abx <- as.character(test_df_abx$sbi_present)

# Candidate thresholds
thresholds <- seq(0.01, 0.15, by = 0.001)

# Compute NPV at each threshold
npv_by_threshold_abx <- purrr::map_dfr(thresholds, function(ocp_rf) {

  # Same classification logic as your no-abx example:
  # predict "yes" if predicted probability >= threshold, else "no"
  pred_class <- ifelse(pred_prob_abx >= ocp_rf, "yes", "no")

  tn <- sum(pred_class == "no" & truth_abx == "no", na.rm = TRUE)
  fn <- sum(pred_class == "no" & truth_abx == "yes", na.rm = TRUE)
  n_pred_negative <- tn + fn

  tibble::tibble(
    threshold = ocp_rf,
    npv = ifelse(n_pred_negative > 0, tn / n_pred_negative, NA_real_),
    n_pred_negative = n_pred_negative,
    n_true_negative = tn,
    n_false_negative = fn
  )
})

View(npv_by_threshold_abx)








# Convert prob -> class using OOF-selected threshold on P(yes)
rf_pred_class_thr <- ifelse(rf_pred_prob_abx[["yes"]] <= threshold_yes_abx, "no", "yes")
rf_pred_class_thr <- factor(rf_pred_class_thr, levels = levels(test_df_abx$sbi_present))

rf_cm_yes_abx_test <- confusionMatrix(
  data      = rf_pred_class_thr,
  reference = test_df_abx$sbi_present,
  mode      = "everything",
  positive  = "yes"
)

npv_retro_yes_abx_test <- unname(rf_cm_yes_abx_test$byClass["Neg Pred Value"])




retro_test_auroc_yes_abx <- as.numeric(
  pROC::auc(
    pROC::roc(
      response  = test_df_abx$sbi_present,
      predictor = rf_pred_prob_abx[["yes"]],
      levels    = c("no", "yes"),
      direction = "<",
      quiet     = TRUE
    )
  )
)

# Store AUPRC in test set also
# make a data frame with truth + positive-class prob
test_pr_df_abx <- tibble(
  truth = test_df_abx$sbi_present,
  p_yes = rf_pred_prob_abx[, "yes"]   # <-- positive class prob (adjust if your level name differs)
) %>%
  mutate(truth = factor(truth, levels = c("yes","no")))

retro_test_auprc_yes_abx <- yardstick::pr_auc(test_pr_df_abx, truth = truth, p_yes) %>% # 0.760 in test set, prevalence 0.48
  pull(.estimate)

# Store the training set AURPC
# AUPRC using PRROC (class0 = positive class scores)
oof_abx <- rf_model_abx$pred %>% mutate(obs = factor(obs, levels = c("yes","no")))
retro_train_auprc_yes_abx <- yardstick::pr_auc(
  oof_abx,
  truth = obs,
  yes        # <-- bare column name, no `estimate =`
) %>%
  pull(.estimate) #0.757 AUPRC, prevalence in training set of SBI is 0.505




# Store OOF probsd
oof_preds_abx <- rf_model_abx$pred %>%
  dplyr::select(rowIndex, obs, yes) %>%
  dplyr::rename(
    sbi_present = obs,
    oof_prob = yes
  )

train_df_with_oof_abx <- train_df_abx %>%
  dplyr::mutate(rowIndex = dplyr::row_number()) %>%
  dplyr::left_join(oof_preds_abx, by = "rowIndex")


# Fix and store train_df_abx_with_oof
train_df_with_oof_abx <- train_df_with_oof_abx %>% dplyr::select(-sbi_present.y) %>% rename(sbi_present = sbi_present.x)
train_df_with_oof_abx <- train_df_with_oof_abx %>% rename(model_score = oof_prob) %>% dplyr::select(-rowIndex)

# write_csv(x = train_df_with_oof_abx, file = "/home/martinbl/sbi_blake/validation_files/rf_ae_40_training_updated_2_4_26.csv")


test_preds_abx <- tibble(
  rowIndex    = seq_len(nrow(test_df_abx)),
  sbi_present = test_df_abx$sbi_present,     # yes/no factor
  model_score = rf_pred_prob_abx[, "yes"]# P(SBI = yes)
)

roc_obj_abx <- roc(
  response  = test_preds_abx$sbi_present,
  predictor = test_preds_abx$model_score,
  levels = c("no", "yes"),   # sets "yes" as the positive class
  direction = "<"
)

auc(roc_obj_abx) # AUROC 0.77

# Now calculate NPV at OOF-selected threshold

test_preds_npv_abx <- test_preds_abx %>%
  mutate(
    pred_class = if_else(model_score <= threshold_yes_abx, "no", "yes"),
    pred_class = factor(pred_class, levels = c("no", "yes")),
    sbi_present = factor(sbi_present, levels = c("no", "yes"))
  )

# NPV = TN / (TN + FN)
npv_retro_yes_abx <- with(test_preds_npv_abx,
                          sum(pred_class == "no" & sbi_present == "no") /
                            sum(pred_class == "no"))

npv_retro_yes_abx # NPV 0.83


test_df_with_pred_abx <- test_df_abx %>%
  mutate(rowIndex = dplyr::row_number()) %>%
  left_join(test_preds_abx %>% select(rowIndex, model_score), by = "rowIndex") %>%
  select(-rowIndex) %>%
  relocate(case_id, sbi_present, model_score)

# Sanity checks
stopifnot(nrow(test_df_with_pred_abx) == nrow(test_df_abx))
stopifnot(!anyNA(test_df_with_pred_abx$model_score))

rf_df_abx <- test_df_with_pred_abx
