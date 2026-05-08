setwd("~/Projects/gentropy/notebooks/gentropy_paper/data/")

library(data.table)
library(arrow)
library(xgboost)

fm=open_dataset("2506_fm_matrix/")
fm=data.frame(fm)
fm$id=paste(fm$studyLocusId,fm$geneId,sep="_")

feature_cols=colnames(fm)
feature_cols=feature_cols[-which(feature_cols%in%c("geneId","studyLocusId","isProteinCoding","id"))]

fm_matrix <- as.matrix(fm[, feature_cols])
dfm <- xgb.DMatrix(data = fm_matrix)

##### First set

training_set_list=c("250624_training_set",
                    "2506_train_set","20250624_full_EGL_max2_no_inter_training_set")

params <- list(
  objective = "binary:logistic",   
  eval_metric = "aucpr",         
  eta = 0.05,                   
  max_depth = 4,                  
  min_child_weight = 20,         
  subsample = 0.7,                
  colsample_bytree = 0.7,     
  lambda = 20,                 
  alpha = 20,
  scale_pos_weight = 0.3
)


for (tr in training_set_list){
  print(tr)
  training_set=open_dataset(tr)
  training_set=data.frame(training_set)
  
  training_set=training_set[training_set$studyLocusId%in%fm$studyLocusId,]

  training_set$label=0
  training_set$label[training_set$goldStandardSet=="positive"]=1
  
  training_set$id=paste(training_set$studyLocusId,training_set$geneId,sep="_")

  ind=match(training_set$id,fm$id)
  table(fm$id[ind]==training_set$id)
  
  X_train <- as.matrix(fm[ind,feature_cols])
  y_train <- training_set$label
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 300
  )
  
  predictions <- predict(model, dfm)
  
  result=cbind(studyLocusId=fm$studyLocusId,geneId=fm$geneId,score=predictions)
  result=as.data.frame(result)
  result=result[result$score>=0.05,]
  
  fwrite(x=result,file=paste0("./l2g_predictions/",tr,".tsv"))
}  

##### Second set

training_set_list=c("2503_training_set",
                    "20250624_full_EGL_max2_string07_training_set",
                    "20250624_full_EGL_max2_string08_dist_training_set",
                    "20250625_restrict_EGL_v1_max2_string08_dist_training_set",
                    "20250625_2506_095_otg_chembl_dist_string08_max2_v1")

training_set_list="20250625_gentropy_paper_v1"

params <- list(
  objective = "binary:logistic",   
  eval_metric = "aucpr",         
  eta = 0.05,                   
  max_depth = 4,                  
  min_child_weight = 20,         
  subsample = 0.7,                
  colsample_bytree = 0.7,     
  lambda = 20,                 
  alpha = 20,
  scale_pos_weight = 0.7
)


for (tr in training_set_list){
  print(tr)
  training_set=open_dataset(tr)
  training_set=data.frame(training_set)
  
  training_set=training_set[training_set$studyLocusId%in%fm$studyLocusId,]
  
  training_set$label=0
  training_set$label[training_set$goldStandardSet=="positive"]=1
  
  training_set$id=paste(training_set$studyLocusId,training_set$geneId,sep="_")
  
  ind=match(training_set$id,fm$id)
  table(fm$id[ind]==training_set$id)
  
  X_train <- as.matrix(fm[ind,feature_cols])
  y_train <- training_set$label
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 300
  )
  
  predictions <- predict(model, dfm)
  
  result=cbind(studyLocusId=fm$studyLocusId,geneId=fm$geneId,score=predictions)
  result=as.data.frame(result)
  result=result[result$score>=0.05,]
  
  fwrite(x=result,file=paste0("./l2g_predictions/",tr,".tsv"))
}  

##### Third set

training_set_list=c("2503_training_set",
                    "20250624_full_EGL_max2_string07_training_set",
                    "20250624_full_EGL_max2_string08_dist_training_set",
                    "20250625_restrict_EGL_v1_max2_string08_dist_training_set",
                    "20250625_2506_095_otg_chembl_dist_string08_max2_v1")

training_set_list="20250625_gentropy_paper_v1"

params <- list(
  objective = "binary:logistic",   # For binary classification
  eval_metric = "aucpr",         
  eta = 0.05,                   
  max_depth =5,                  
  min_child_weight = 10,         
  subsample = 0.8,                
  colsample_bytree = 0.8,     
  lambda = 1,                 
  alpha = 1,
  scale_pos_weight = 0.8
)


for (tr in training_set_list){
  print(tr)
  training_set=open_dataset(tr)
  training_set=data.frame(training_set)
  
  training_set=training_set[training_set$studyLocusId%in%fm$studyLocusId,]
  
  training_set$label=0
  training_set$label[training_set$goldStandardSet=="positive"]=1
  
  training_set$id=paste(training_set$studyLocusId,training_set$geneId,sep="_")
  
  ind=match(training_set$id,fm$id)
  table(fm$id[ind]==training_set$id)
  
  X_train <- as.matrix(fm[ind,feature_cols])
  y_train <- training_set$label
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 300
  )
  
  predictions <- predict(model, dfm)
  
  result=cbind(studyLocusId=fm$studyLocusId,geneId=fm$geneId,score=predictions)
  result=as.data.frame(result)
  result=result[result$score>=0.05,]
  
  fwrite(x=result,file=paste0("./l2g_predictions/",tr,"_v2.tsv"))
}  




