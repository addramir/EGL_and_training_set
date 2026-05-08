setwd("Desktop/training_set/")
library(data.table)


training_set_list=c("patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_studies.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_ALL_no_filter.tsv",
                    "patched_training_2503-testrun-1_v5_1_no_exposion.tsv",
                    "patched_training_2503-testrun-1_old_gs.tsv",
                    "training_2503-testrun-1_all_string_liberal_egl.tsv")
test_set="test_patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv"


#### 
#### 
#### 
library(data.table)
library(xgboost)
library(PRROC)
library(pROC)

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
###
### test 1
###

test=fread(test_set,data.table=F)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
y_test=test$GSP
test=test[,feat_colnames]
test_matrix <- as.matrix(test)
dtest <- xgb.DMatrix(data = test_matrix)

#### 
#### model 1
#### 
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

ts=training_set_list[1]

for (ts in c(training_set_list)){
  fm=fread(ts,data.table=F)
  
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  
  ind=which(fm$studyLocusId%in%fm1_cs) 
  train=fm[ind,feat_colnames]
  y_train=fm[ind,"GSP"]
  train_matrix <- as.matrix(train)
  dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
  
  print(ts)
  
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,
    watchlist = list(train = dtrain),
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose=0
  )
  
  probabilities <- predict(xgb_model, dtest)
  pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  
  auc_pr <- pr$auc.integral  # AUC-PR score
  print(paste("AUC-PR:", auc_pr))
  
  auc_result <- roc(y_test, probabilities)
  print(auc_result)
  cat("\n")
  cat("\n")
}



#### 
#### model 2
#### 

params <- list(
  objective = "binary:logistic",
  booster = "gbtree",
  eta = 0.05,
  max_depth = 10,
  min_child_weight = 10,
  subsample = 1,
  colsample_bytree = 1,
  lambda = 1,
  alpha = 0,
  eval_metric = "auc",
  gamma=0
)

for (ts in c(training_set_list)){
  fm=fread(ts,data.table=F)
  
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  
  ind=which(fm$studyLocusId%in%fm1_cs) 
  train=fm[ind,feat_colnames]
  y_train=fm[ind,"GSP"]
  train_matrix <- as.matrix(train)
  dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
  
  print(ts)
 
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,
    watchlist = list(train = dtrain),
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose=0
  )
  
  probabilities <- predict(xgb_model, dtest)
  pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  
  auc_pr <- pr$auc.integral  # AUC-PR score
  print(paste("AUC-PR:", auc_pr))
  
  auc_result <- roc(y_test, probabilities)
  print(auc_result)
}




###
### test 2
###
test_set="test_patched_training_2503-testrun-1_old_gs.tsv"
test=fread(test_set,data.table=F)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
y_test=test$GSP
test=test[,feat_colnames]
test_matrix <- as.matrix(test)
dtest <- xgb.DMatrix(data = test_matrix, label = y_test)

#### 
#### model 2
#### 
params <- list(
  objective = "binary:logistic",
  booster = "gbtree",
  eta = 0.05,
  max_depth = 10,
  min_child_weight = 10,
  subsample = 1,
  colsample_bytree = 1,
  lambda = 1,
  alpha = 0,
  eval_metric = "auc",
  gamma=0
)

for (ts in c(training_set_list)){
  fm=fread(ts,data.table=F)
  
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  
  ind=which(fm$studyLocusId%in%fm1_cs) 
  train=fm[ind,feat_colnames]
  y_train=fm[ind,"GSP"]
  train_matrix <- as.matrix(train)
  dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
  
  print(ts)
  
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,
    watchlist = list(train = dtrain),
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose=0
  )
  
  probabilities <- predict(xgb_model, dtest)
  pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  
  auc_pr <- pr$auc.integral  # AUC-PR score
  print(paste("AUC-PR:", auc_pr))
  
  auc_result <- roc(y_test, probabilities)
  print(auc_result)
}



###
### test 3
###
test_set="test_patched_training_2503-testrun-1_v5_1_no_exposion.tsv"
test=fread(test_set,data.table=F)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
y_test=test$GSP
test=test[,feat_colnames]
test_matrix <- as.matrix(test)
dtest <- xgb.DMatrix(data = test_matrix, label = y_test)

#### 
#### model 1
#### 
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

for (ts in c(training_set_list)){
  fm=fread(ts,data.table=F)
  
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  
  ind=which(fm$studyLocusId%in%fm1_cs) 
  train=fm[ind,feat_colnames]
  y_train=fm[ind,"GSP"]
  train_matrix <- as.matrix(train)
  dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
  
  print(ts)
  
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,
    watchlist = list(train = dtrain),
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose=0
  )
  
  probabilities <- predict(xgb_model, dtest)
  pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  
  auc_pr <- pr$auc.integral  # AUC-PR score
  print(paste("AUC-PR:", auc_pr))
  
  auc_result <- roc(y_test, probabilities)
  print(auc_result)
}
