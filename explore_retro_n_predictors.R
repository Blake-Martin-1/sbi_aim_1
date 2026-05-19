### Code to explore how different number of inputs affect the AUROC when applied to the retrospective test set ###

# First divide the rf_df into the train and test data frames by random 80/20 split
set.seed(seed = 32313)

train_size <- floor(nrow(rf_df) * 0.8)

train_index <- sample.int(
  n = nrow(rf_df),
  size = train_size,
  replace = FALSE
)

train_df <- rf_df[train_index, ]
test_df  <- rf_df[-train_index, ]

predictors_full_retro <- names(
  train_df %>%
    dplyr::select(-case_id, -sbi_present)
)

rf_train_df <- train_df %>%
  dplyr::select(sbi_present, dplyr::all_of(predictors_full_retro))

rf_model_trial <- caret::train(
  sbi_present ~ .,
  data = rf_train_df,
  method = "ranger",
  tuneLength = 10,
  na.action = na.omit,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  trControl = caret::trainControl(
    method = "cv",
    number = 5,
    summaryFunction = caret::twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE
  )
)