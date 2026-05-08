setwd("Desktop/training_set/")
library(data.table)

#x=fread("training_2503-testrun-1_all_string_liberal_studies.tsv",data.table=F)
#x=cbind(x,sid_gid=paste(x$studyLocusId,x$geneId,sep="_"))

x=fread("training_2503-testrun-1_all_string_liberal_variants.tsv",data.table=F)

#x1=cbind(x1,sid_gid=paste(x1$studyLocusId,x1$geneId,sep="_"))
#table(x1$sid_gid%in%x$sid_gid)

feat_colnames=c("credibleSetConfidence","distanceFootprintMean",
"distanceFootprintMeanNeighbourhood","distanceSentinelFootprint",
"distanceSentinelFootprintNeighbourhood","distanceSentinelTss",
"distanceSentinelTssNeighbourhood","distanceTssMean",
"distanceTssMeanNeighbourhood","eQtlColocClppMaximum",
"eQtlColocClppMaximumNeighbourhood","eQtlColocH4Maximum",
"eQtlColocH4MaximumNeighbourhood","geneCount500kb",
"isProteinCoding","pQtlColocClppMaximum",
"pQtlColocClppMaximumNeighbourhood","pQtlColocH4Maximum",
"pQtlColocH4MaximumNeighbourhood","proteinGeneCount500kb",
"sQtlColocClppMaximum","sQtlColocClppMaximumNeighbourhood",
"sQtlColocH4Maximum","sQtlColocH4MaximumNeighbourhood",
"vepMaximum","vepMaximumNeighbourhood",
"vepMean","vepMeanNeighbourhood")

fm=x[,feat_colnames]
y=x$GSP

gene_efo=unique(x$geneid_efo[x$GSP==1])


#WAY 1
set.seed(43)
prop=0.2
#ind_train=sample(1:length(y),length(y)*(1-prop))
ind_efo=sample(1:length(gene_efo),length(gene_efo)*(1-prop))
cs_in_train=unique(x$studyLocusId[x$geneid_efo%in%gene_efo[ind_efo]])
ind_train=which(x$studyLocusId%in%cs_in_train)

train=fm[ind_train,]
y_train=y[ind_train]

test=fm[-ind_train,]
y_test=y[-ind_train]


####
library(xgboost)
library(PRROC)
library(pROC)

# Count classes
num_negative <- sum(y_train == 0)
num_positive <- sum(y_train == 1)

# Compute scale_pos_weight
scale_weight <- num_negative / num_positive
print(paste("Scale Pos Weight:", scale_weight))

train_matrix <- as.matrix(train)
test_matrix <- as.matrix(test)

dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
dtest <- xgb.DMatrix(data = test_matrix, label = y_test)

params <- list(
  objective = "binary:logistic",  # Binary classification
  booster = "gbtree",
  eta = 0.01,                      # Learning rate
  max_depth = 5,                  # Tree depth
  min_child_weight = 10,            # Min data per leaf
  subsample = 0.8,                 # Sample ratio of training data
  colsample_bytree = 1,          # Feature sampling
  gamma = 0.1,                       # Min loss reduction
  eval_metric = "aucpr",           # AUC-PR for imbalanced datasets
  lambda=5
)

# Train model with early stopping
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  watchlist = list(train = dtrain,test=dtest),
  early_stopping_rounds = 10,
  maximize = TRUE,
  alpha=0
)

# Predict probabilities
probabilities <- predict(xgb_model, dtest)

# Compute AUC-PR (Recommended for Imbalanced Data)

pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)

auc_pr <- pr$auc.integral  # AUC-PR score
print(paste("AUC-PR:", auc_pr))

# Compute ROC-AUC (Optional)

auc_result <- roc(y_test, probabilities)
print(auc_result)


#######
param_grid <- expand.grid(
  eta = c(0.005, 0.01, 0.05),           # Learning rates to test
  max_depth = c(3, 5, 10),               # Tree depths
  min_child_weight = c(5, 10, 20),      # Min observations per leaf
  subsample = c(0.8,1),                # Sample fraction of data
  lambda = c(1, 5, 10),                 # L2 Regularization
  alpha = c(0, 0.5, 1),                  # L1 Regularization
  gamma=c(0,0.1,1,5)
)

best_model <- NULL
best_aucpr <- NULL

# Loop through each parameter combination
for (i in 1:nrow(param_grid)) {
  
  # Extract parameters for current iteration
  params <- list(
    objective = "binary:logistic",
    booster = "gbtree",
    eta = param_grid$eta[i],
    max_depth = param_grid$max_depth[i],
    min_child_weight = param_grid$min_child_weight[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = 1,
    lambda = param_grid$lambda[i],
    alpha = param_grid$alpha[i],
    eval_metric = "auc",
    gamma=param_grid$gamma[i]
  )
  
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,
    watchlist = list(train = dtrain,test=dtest),
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose=0
  )
  
  # Predict probabilities
  #probabilities <- predict(xgb_model, dtest)
  #pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  #auc_pr <- pr$auc.integral  # AUC-PR score
  #print(paste(i,"AUC-PR:", auc_pr))
  #best_aucpr=c(best_aucpr,auc_pr)
  
  #print(paste(i,"AUC-PR:", max(xgb_model$evaluation_log$test_aucpr)))
  #best_aucpr=c(best_aucpr,max(xgb_model$evaluation_log$test_aucpr))
  print(paste(i,"AUC:", max(xgb_model$evaluation_log$test_auc)))
  best_aucpr=c(best_aucpr,max(xgb_model$evaluation_log$test_auc))

}

summary(best_aucpr)
i=which.max(best_aucpr)

summary(lm(best_aucpr~as.matrix(param_grid)))


params <- list(
  objective = "binary:logistic",
  booster = "gbtree",
  eta = 0.05,
  max_depth = 5,
  min_child_weight = 10,
  subsample = 0.8,
  colsample_bytree = 1,
  lambda = 10,
  alpha = 0,
  eval_metric = "aucpr",
  gamma=0
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  watchlist = list(train = dtrain),
  early_stopping_rounds = 10,
  maximize = TRUE,
  verbose=0
)
auc_pr <- pr$auc.integral  # AUC-PR score
print(paste("AUC-PR:", auc_pr))

# Compute ROC-AUC (Optional)

auc_result <- roc(y_test, probabilities)
print(auc_result)

ap_result <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
print(paste("Average Precision (AP):", ap_result$auc.integral))

#WAY 2
cs_test=fread("study_loci_to_test.txt",data.table=F,header=F)
cs_test=cs_test[,1]
ind_test=which(x$studyLocusId%in%cs_test)


train=fm[-ind_test,]
y_train=y[-ind_test]

test=fm[ind_test,]
y_test=y[ind_test]

####
library(brglm2)
library(pROC)

firth_model <- glm(y_train ~ ., data = train, family = binomial(link = "logit"), method = brglmFit)

probabilities <- predict(firth_model, newdata = test, type = "response")
auc_result <- roc(y_test, probabilities)
print(auc_result)


library(nnet)
nn_model <- nnet(y_train ~ ., data = train, size = 10, linout = FALSE, entropy = TRUE, maxit = 500)
probabilities <- predict(nn_model, newdata = test, type = "raw")
auc_result <- roc(y_test, probabilities)
print(auc_result)

library(neuralnet)
nn_model <- neuralnet(y_train ~ ., data = train, hidden = c(10, 5), linear.output = FALSE)
nn_results <- compute(nn_model, test)
probabilities <- nn_results$net.result
auc_result <- roc(y_test, probabilities)
print(auc_result)



######
library(xgboost)

# Convert data to matrix
train_matrix <- as.matrix(train)
test_matrix <- as.matrix(test)

# Create DMatrix objects
dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
dtest <- xgb.DMatrix(data = test_matrix, label = y_test)

# Set parameters to match Python settings
params <- list(
  objective = "binary:logistic",   # Equivalent to log_loss
  booster = "gbtree",              # Tree-based model
  eta = 0.1,                       # Learning rate
  max_depth = 10,                  # Tree depth
  min_child_weight = 5,            # Equivalent to min_samples_leaf
  gamma = 0,                       # Equivalent to min_impurity_decrease
  subsample = 1,                   # Use all data in each boosting round
  colsample_bytree = 1,            # Equivalent to max_features=null (all features used)
  nthread = parallel::detectCores() # Use all available cores
)

# Train XGBoost Model
xgb_model <- xgboost(
  data = dtrain,
  params = params,
  nrounds = 500,       # Equivalent to n_estimators
  early_stopping_rounds = NULL,  # No equivalent to n_iter_no_change in R
  verbose = 0           # Silent mode (Python’s verbose=0)
)

# Predict probabilities
probabilities <- predict(xgb_model, dtest)

# Compute AUC
library(pROC)
library(PRROC)
auc_result <- roc(y_test, probabilities)
print(auc_result)
mAP_result <- mAP(probabilities, y_test)
print(paste("mAP Score:", mAP_result))

param_grid <- expand.grid(
  max_depth = c(10,20),
  eta = c(0.01, 0.05, 0.1, 0.2, 0.3),
  min_child_weight = c(3),
  subsample = c(0.8, 1),
  colsample_bytree = c(0.8, 1),
  gamma = c(0, 0.1)
)

best_auc <- NULL
best_params <- NULL

for (i in 1:nrow(param_grid)) {
  print(i)
  params <- list(
    objective = "binary:logistic",
    booster = "gbtree",
    eta = param_grid$eta[i],
    max_depth = param_grid$max_depth[i],
    min_child_weight = param_grid$min_child_weight[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma = param_grid$gamma[i],
    eval_metric = "auc"
  )
  
  
  xgb_model <- xgboost(
    data = dtrain,
    params = params,
    nrounds = 100,       # Equivalent to n_estimators
    early_stopping_rounds = NULL,  # No equivalent to n_iter_no_change in R
    verbose = 0           # Silent mode (Python’s verbose=0)
  )

  probabilities <- predict(xgb_model, dtest)
  
  auc_result <- roc(y_test, probabilities)
  
  best_auc=c(best_auc,auc_result$auc[1])
}

print(best_params)
print(best_auc)


