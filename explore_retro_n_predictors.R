### Code to explore how different number of inputs affect the AUROC when applied to the retrospective test set ###

# First divide the rf_df into the train and test data frames by random 80/20 split
set.seed(seed = 32313)
train_size <- floor(nrow(rf_df) * 0.8)
sample <- sample.int(n = nrow(rf_df), size = train_size, replace = FALSE)
train_df <- rf_df[sample, ] #8228 pts
test_df <- rf_df[-sample, ] # 2058 patients


# # Now will train / creat the RF model for abx-unexposed patients
# Train an initial RF model across all mtry and gini type values
rf_model_trial <- train(sbi_present ~ ., tuneLength = 10, data = train_df, method = "ranger", na.action = na.omit, importance = 'impurity',
                        respect.unordered.factors = TRUE,
                        trControl = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE,
                                                 verboseIter = TRUE))