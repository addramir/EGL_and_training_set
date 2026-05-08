setwd("~/Projects/gentropy/notebooks/gentropy_paper/data/")


###
fun_for_xgboost=function(traing_data_set,parameters,nrounds=500){
  # 5. Train the model
  model <- xgb.train(
    params = parameters,
    data = traing_data_set,
    nrounds = nrounds
  )
  
  predictions <- predict(model, dfm)
  ind=which(predictions>=0.5)
  result=cbind(studyLocusId=fm$studyLocusId[ind],geneId=fm$geneId[ind],l2g=predictions[ind])
  result=as.data.frame(result)
  print(paste("total number of studyLocusId with l2g>0.5:",length(unique(result$studyLocusId))))
  lll=table(result$studyLocusId)
  print(paste("total number of studyLocusId with one gene with l2g>0.5:",sum(lll==1)))
  print(paste("total number of studyLocusId with two genes with l2g>0.5:",sum(lll==2)))
  print(paste("total number of studyLocusId with more than two genes with l2g>0.5:",sum(lll>2)))
  
  predictions_2503 <- predict(model, dfm_2503)
  ind=which(predictions_2503>=0.5)
  result=cbind(studyLocusId=fm_2503$studyLocusId[ind],geneId=fm_2503$geneId[ind],l2g=predictions_2503[ind])
  result=as.data.frame(result)
  print(paste("total number of studyLocusId on 2503 fm with l2g>0.5:",length(unique(result$studyLocusId))))
  lll=table(result$studyLocusId)
  print(paste("total number of studyLocusId on 2503 fm with one gene with l2g>0.5:",sum(lll==1)))
  print(paste("total number of studyLocusId on 2503 fm with two genes with l2g>0.5:",sum(lll==2)))
  print(paste("total number of studyLocusId on 2503 fm with more than two genes with l2g>0.5:",sum(lll>2)))
  
  ind=which(l2g_2503$id%in%fm_2503$id)
  l2g_tmp=l2g_2503[ind,]
  ind=match(l2g_tmp$id,fm_2503$id)
  
  print(paste("cor with l2g 2503:",cor(l2g_tmp$score,predictions_2503[ind])))
  
  ind1=sample(1:length(l2g_tmp$score),10000)
  plot(l2g_tmp$score[ind1],predictions[ind][ind1])
  

  print(paste("min predictions_2503:",min(predictions_2503),"max:", max(predictions_2503)))
  return(predictions_2503)
}
###


library(data.table)
library(arrow)
library(xgboost)

fm=open_dataset("2506_fm_matrix/")
fm=data.frame(fm)

fm_2503=open_dataset("2503_fm_matrix/")
fm_2503=data.frame(fm_2503)
fm_2503$id=paste(fm_2503$studyLocusId,fm_2503$geneId,sep="_")

l2g_2503=open_dataset("2503_l2g/")
l2g_2503=data.frame(l2g_2503)
l2g_2503=l2g_2503[,c("studyLocusId","geneId","score")]
l2g_2503$id=paste(l2g_2503$studyLocusId,l2g_2503$geneId,sep="_")

training_set=open_dataset("2503_training_set/")
training_set=data.frame(training_set)

training_set=open_dataset("250624_training_set/")
training_set=data.frame(training_set)

training_set=open_dataset("2506_train_set/")
training_set=data.frame(training_set)

training_set=open_dataset("20250624_full_EGL_max2_no_inter_training_set/")
training_set=data.frame(training_set)

training_set=open_dataset("20250624_full_EGL_max2_string07_training_set/")
training_set=data.frame(training_set)

training_set=open_dataset("20250624_full_EGL_max2_string08_dist_training_set/")
training_set=data.frame(training_set)


table(training_set$studyLocusId%in%fm$studyLocusId)
#FALSE   TRUE 
#2604 100427 

training_set=training_set[training_set$studyLocusId%in%fm$studyLocusId,]
table(training_set$goldStandardSet)
#negative positive 
#92076     8351 

training_set$label=0
training_set$label[training_set$goldStandardSet=="positive"]=1

training_set$id=paste(training_set$studyLocusId,training_set$geneId,sep="_")
fm$id=paste(fm$studyLocusId,fm$geneId,sep="_")
table(training_set$id%in%fm$id)

feature_cols=colnames(fm)
feature_cols=feature_cols[-which(feature_cols%in%c("geneId","studyLocusId","isProteinCoding","id"))]

table(fm$pQtlColocClppMaximum>=0.01 & 
        (fm$distanceSentinelFootprintNeighbourhood==1|fm$distanceSentinelTssNeighbourhood==1))
#FALSE     TRUE 
#10589844    33527
####

fm_mod=fm[fm$id%in%training_set$id,]
fm_mod[fm_mod[,"vepMaximum"]<=0.33,"vepMaximum"]=0
x=c("eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum","vepMaximum")
ll=apply(fm_mod[,x],MAR=1,sum)
fm_mod=cbind(fm_mod,sum_col=ll)
ind=match(training_set$id,fm_mod$id)
table(training_set$id==fm_mod$id[ind])
fm_mod=fm_mod[ind,]
table(fm_mod$id==training_set$id)
fm_mod=cbind(fm_mod,gsp=training_set$label)
poss=fm_mod[fm_mod$gsp==1,]
table(poss$sum_col==0)
table(poss$distanceSentinelFootprint>=0.1)

table(poss$sum_col==0 & 
        !(poss$distanceSentinelFootprintNeighbourhood==1|poss$distanceSentinelTssNeighbourhood==1))

ind=which(poss$sum_col==0 & 
            !(poss$distanceSentinelFootprintNeighbourhood==1|poss$distanceSentinelTssNeighbourhood==1))
studyLoci_to_exclude=poss$studyLocusId[ind]
training_set2=training_set[!(training_set$studyLocusId%in%studyLoci_to_exclude),]

ind=which(poss$sum_col>0.01 & 
            (poss$distanceSentinelFootprintNeighbourhood==1|poss$distanceSentinelTssNeighbourhood==1))
studyLoci_to_include=poss$studyLocusId[ind]
training_set3=training_set[training_set$studyLocusId%in%studyLoci_to_include,]

ind=which(poss$sum_col==0 & 
            !(poss$distanceSentinelFootprintNeighbourhood==1|poss$distanceSentinelTssNeighbourhood==1))
studyLoci_to_include=poss$studyLocusId[ind]
training_set4=training_set[training_set$studyLocusId%in%studyLoci_to_include,]

ind=which(poss$pQtlColocClppMaximum>=0.01 & 
                    (poss$distanceSentinelFootprintNeighbourhood==1|poss$distanceSentinelTssNeighbourhood==1))
studyLoci_to_include=poss$studyLocusId[ind]
training_set5=training_set[training_set$studyLocusId%in%studyLoci_to_include,]


unique_genes=unique(poss$geneId)
genes_to_select=sample(unique_genes,size = length(unique_genes)*0.8)
studyLoci_to_include=unique(poss$studyLocusId[poss$geneId%in%genes_to_select])
training_set7=training_set[training_set$studyLocusId%in%studyLoci_to_include,]

ind=which(training_set$goldStandardSet=="negative")
ind=sample(ind,length(ind)*0.5)
training_set8=training_set[-ind,]

############
ind=match(training_set$id,fm$id)
table(fm$id[ind]==training_set$id)

X_train <- as.matrix(fm[ind,feature_cols])
y_train <- training_set$label

############
ind=match(training_set2$id,fm$id)
table(fm$id[ind]==training_set2$id)

X_train2 <- as.matrix(fm[ind,feature_cols])
y_train2 <- training_set2$label
############
ind=match(training_set3$id,fm$id)
table(fm$id[ind]==training_set3$id)

X_train3 <- as.matrix(fm[ind,feature_cols])
y_train3 <- training_set3$label
############
ind=match(training_set4$id,fm$id)
table(fm$id[ind]==training_set4$id)

X_train4 <- as.matrix(fm[ind,feature_cols])
y_train4 <- training_set4$label
############
ind=match(training_set5$id,fm$id)
table(fm$id[ind]==training_set5$id)

X_train5 <- as.matrix(fm[ind,feature_cols])
y_train5 <- training_set5$label
############
ind=match(training_set$id,fm$id)
table(fm$id[ind]==training_set$id)

X_train6 <- as.matrix(fm[ind,feature_cols])
y_train6 <- sample(training_set$label)
##########
ind=match(training_set7$id,fm$id)
table(fm$id[ind]==training_set7$id)

X_train7 <- as.matrix(fm[ind,feature_cols])
y_train7 <- training_set7$label
############
ind=match(training_set8$id,fm$id)
table(fm$id[ind]==training_set8$id)

X_train8 <- as.matrix(fm[ind,feature_cols])
y_train8 <- training_set8$label
############

# Convert to DMatrix
dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtrain2 <- xgb.DMatrix(data = X_train2, label = y_train2)
dtrain3 <- xgb.DMatrix(data = X_train3, label = y_train3)
dtrain4 <- xgb.DMatrix(data = X_train4, label = y_train4)
dtrain5 <- xgb.DMatrix(data = X_train5, label = y_train5)
dtrain6 <- xgb.DMatrix(data = X_train6, label = y_train6)
dtrain7 <- xgb.DMatrix(data = X_train7, label = y_train7)
dtrain8 <- xgb.DMatrix(data = X_train8, label = y_train8)

fm_matrix <- as.matrix(fm[, feature_cols])  # assumes fm has same features
dfm <- xgb.DMatrix(data = fm_matrix)

fm_matrix_2503 <- as.matrix(fm_2503[, feature_cols])  # assumes fm has same features
dfm_2503 <- xgb.DMatrix(data = fm_matrix_2503)


# 4. Set parameters
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




l=fun_for_xgboost(traing_data_set = dtrain,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain2,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain3,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain4,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain5,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain6,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain7,parameters = params,nrounds = 300)
l=fun_for_xgboost(traing_data_set = dtrain8,parameters = params,nrounds = 300)





l2g_2503_filtered=l2g_2503[l2g_2503$score>=0.5,]
print(paste("total number of studyLocusId with l2g>0.5:",length(unique(l2g_2503_filtered$studyLocusId))))
lll=table(l2g_2503_filtered$studyLocusId)
print(paste("total number of studyLocusId with one gene with l2g>0.5:",sum(lll==1)))
print(paste("total number of studyLocusId with two genes with l2g>0.5:",sum(lll==2)))
print(paste("total number of studyLocusId with more than two genes with l2g>0.5:",sum(lll>2)))

# Build the data frame
ind=which(l2g_2503$id%in%fm_2503$id)
l2g_tmp=l2g_2503[ind,]
ind=match(l2g_tmp$id,fm_2503$id)
df <- data.frame(
  l2g_2503 = l2g_tmp$score,
  new = l[ind]
)
summary(lm(df$new~df$l2g_2503))


