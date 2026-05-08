setwd("Desktop/training_set/")
library(data.table)

x=fread("v5_1_no_explosion_fm.csv/part-00000-e57891c7-ed63-4ece-add4-0c0445dcfae8-c000.csv",data.table=F)

dim(x)
x=na.omit(x)
dim(x)
x$GSP=0
x$GSP[x$goldStandardSet=="positive"]=1
table(x$GSP)
x$geneid_efo=paste(x$geneId,x$efo_terms,sep="_")

write.table(x,file="training_2503-testrun-1_v5_1_no_exposion.tsv",sep="\t",quote=F,col.names = T,row.names = F)



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

prop=0.2
#ind_train=sample(1:length(y),length(y)*(1-prop))
ind_efo=sample(1:length(gene_efo),length(gene_efo)*(1-prop))
cs_in_train=unique(x$studyLocusId[x$geneid_efo%in%gene_efo[ind_efo]])
ind_train=which(x$studyLocusId%in%cs_in_train)


train=fm[ind_train,]
y_train=y[ind_train]

test=fm[-ind_train,]
y_test=y[-ind_train]

study_loci_to_test=unique(x$studyLocusId[-ind_train])
write.table(study_loci_to_test,file="study_loci_to_test.txt",sep="\t",quote=FALSE,row.names=F,col.names=F)

#### xgbost
library(pROC)

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
  nrounds = 100,       # Equivalent to n_estimators
  early_stopping_rounds = NULL,  # No equivalent to n_iter_no_change in R
  verbose = 0           # Silent mode (Python’s verbose=0)
)

# Predict probabilities
probabilities <- predict(xgb_model, dtest)

# Compute AUC
library(pROC)
auc_result <- roc(y_test, probabilities)
print(auc_result)


param_grid <- expand.grid(
  max_depth = c(4, 6, 8, 10),
  eta = c(0.01, 0.05, 0.1, 0.2),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.8, 1),
  colsample_bytree = c(0.8, 1),
  gamma = c(0, 0.1, 0.5)
)

best_auc <- 0
best_params <- NULL

for (i in 1:nrow(param_grid)) {
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
  
  cv_model <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 100,            # Number of boosting iterations
    nfold = 5,                # 5-fold cross-validation
    stratified = TRUE,
    showsd = TRUE,
    early_stopping_rounds = 10,
    maximize = TRUE,
    verbose = 0               # Set to 1 to see progress
  )
  
  mean_auc <- max(cv_model$evaluation_log$test_auc_mean)
  
  if (mean_auc > best_auc) {
    best_auc <- mean_auc
    best_params <- params
  }
}

print(best_params)
print(best_auc)


