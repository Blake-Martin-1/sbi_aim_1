### Evaluate how predictor count affects AUROC for SBI RF model ###

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
  library(doParallel)
  library(pROC)
})

#-----------------------------
# 0) Train / test split
#-----------------------------
set.seed(32313)
train_size <- floor(nrow(rf_df) * 0.8)
idx <- sample.int(n = nrow(rf_df), size = train_size, replace = FALSE)
train_df <- rf_df[idx, ]
test_df <- rf_df[-idx, ]

# Ensure outcome is factor and set positive class first for twoClassSummary
if (!is.factor(train_df$sbi_present)) {
  train_df$sbi_present <- factor(train_df$sbi_present)
  test_df$sbi_present <- factor(test_df$sbi_present)
}

guess_positive_class <- function(x) {
  lv <- levels(x)
  hit <- grep("^(yes|y|1|true|positive|pos|sbi)$", lv, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) hit[1] else lv[1]
}

positive_class <- guess_positive_class(train_df$sbi_present)
train_df$sbi_present <- relevel(train_df$sbi_present, ref = positive_class)
test_df$sbi_present <- factor(test_df$sbi_present, levels = levels(train_df$sbi_present))

#-----------------------------
# 1) Parallel backend (n - 1 cores)
#-----------------------------
n_cores <- parallel::detectCores(logical = TRUE)
workers <- max(1, n_cores - 1)
cl <- parallel::makePSOCKcluster(workers)
doParallel::registerDoParallel(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  foreach::registerDoSEQ()
}, add = TRUE)

#-----------------------------
# 2) Repeated CV for tuning on train only
#-----------------------------
ctrl_repeated <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 2,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  allowParallel = TRUE,
  verboseIter = FALSE
)

# Fast but robust random search over RF hyperparameters
set.seed(32313)
rf_tuned <- train(
  sbi_present ~ .,
  data = train_df,
  method = "ranger",
  metric = "ROC",
  trControl = ctrl_repeated,
  tuneLength = 20,
  importance = "permutation",
  respect.unordered.factors = "order",
  num.threads = 1,
  na.action = na.omit
)

best_params <- rf_tuned$bestTune
cat("Best tuned parameters:\n")
print(best_params)

#-----------------------------
# 3) Permutation importance ranking on full train
#-----------------------------
set.seed(32313)
rf_importance <- ranger(
  formula = sbi_present ~ .,
  data = train_df,
  mtry = best_params$mtry,
  splitrule = as.character(best_params$splitrule),
  min.node.size = best_params$min.node.size,
  num.trees = 500,
  probability = TRUE,
  importance = "permutation",
  respect.unordered.factors = "order",
  num.threads = workers,
  na.action = "na.omit"
)

importance_df <- data.frame(
  feature = names(rf_importance$variable.importance),
  importance = as.numeric(rf_importance$variable.importance),
  stringsAsFactors = FALSE
)
importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]
predictors_ranked <- importance_df$feature
p <- length(predictors_ranked)

#-----------------------------
# 4) For k = 1..p, fit RF with top-k features via repeated CV
#-----------------------------
roc_by_k <- numeric(p)

for (k in seq_len(p)) {
  top_k <- predictors_ranked[seq_len(k)]
  train_k <- train_df[, c("sbi_present", top_k), drop = FALSE]

  # mtry cannot exceed number of predictors
  tune_k <- best_params
  tune_k$mtry <- min(tune_k$mtry, k)

  set.seed(32313 + k)
  fit_k <- train(
    sbi_present ~ .,
    data = train_k,
    method = "ranger",
    metric = "ROC",
    trControl = ctrl_repeated,
    tuneGrid = tune_k,
    importance = "none",
    respect.unordered.factors = "order",
    num.threads = 1,
    na.action = na.omit
  )

  roc_by_k[k] <- max(fit_k$results$ROC, na.rm = TRUE)
  cat(sprintf("k=%d/%d, CV AUROC=%.4f\n", k, p, roc_by_k[k]))
}

performance_by_k <- data.frame(
  k = seq_len(p),
  cv_auroc = roc_by_k
)

# Parsimony rule: smallest k within 10% of max AUROC
max_auc <- max(performance_by_k$cv_auroc, na.rm = TRUE)
threshold_auc <- 0.9 * max_auc
selected_k <- min(performance_by_k$k[performance_by_k$cv_auroc >= threshold_auc])
selected_features <- predictors_ranked[seq_len(selected_k)]

cat(sprintf("Max CV AUROC: %.4f\n", max_auc))
cat(sprintf("Parsimony threshold (90%% of max): %.4f\n", threshold_auc))
cat(sprintf("Selected k: %d\n", selected_k))

#-----------------------------
# 5) Refit final model on full train with selected k
#-----------------------------
final_train <- train_df[, c("sbi_present", selected_features), drop = FALSE]
final_tune <- best_params
final_tune$mtry <- min(final_tune$mtry, selected_k)

set.seed(32313)
final_rf <- train(
  sbi_present ~ .,
  data = final_train,
  method = "ranger",
  metric = "ROC",
  trControl = trainControl(method = "none", classProbs = TRUE),
  tuneGrid = final_tune,
  importance = "none",
  respect.unordered.factors = "order",
  num.threads = workers,
  na.action = na.omit
)

#-----------------------------
# 6) Evaluate once on test_df (AUROC + CI)
#-----------------------------
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

cat(sprintf("Final test AUROC: %.4f\n", auc_val))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", auc_ci[1], auc_ci[3]))

# Optional objects for downstream use
results <- list(
  tuned_model = rf_tuned,
  importance = importance_df,
  performance_by_k = performance_by_k,
  selected_k = selected_k,
  selected_features = selected_features,
  final_model = final_rf,
  test_auc = auc_val,
  test_auc_ci = auc_ci
)
