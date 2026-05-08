library(data.table)
library(xgboost)
library(PRROC)
library(pROC)

setwd("Desktop/training_set/")
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

training_set_list=c("patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_studies.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_ALL_no_filter.tsv",
                    "patched_training_2503-testrun-1_v5_1_no_exposion.tsv",
                    "patched_training_2503-testrun-1_old_gs.tsv",
                    "training_2503-testrun-1_all_string_liberal_egl.tsv")

test_set="test_patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv"

pred=fread("20250303_fm_all_rows.csv/part-00000-1974d9aa-c0c2-461b-97d1-592b32a14fea-c000.csv",data.table = F)
pred=pred[pred$isProteinCoding==1,]
length(unique(pred$studyLocusId))
pred_fm=as.matrix(pred[,feat_colnames])
dpred <- xgb.DMatrix(data = pred_fm)


priv_l2g=fread("24_12_l2g.csv/part-00000-709fe0b1-eeba-403d-a226-37a6fb570a43-c000.csv",data.table=FALSE)

test=fread(test_set,data.table=F)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
y_test=test$GSP
test=test[,feat_colnames]
test_matrix <- as.matrix(test)
dtest <- xgb.DMatrix(data = test_matrix)


ts=training_set_list[1]
for (ts in c(training_set_list)){
  fm=fread(ts,data.table=FALSE)
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  
  ind=which(fm$studyLocusId%in%fm1_cs) 
  train=fm[ind,feat_colnames]
  y_train=fm[ind,"GSP"]
  train_matrix <- as.matrix(train)
  dtrain <- xgb.DMatrix(data = train_matrix, label = y_train)
  
  print(ts)

  params <- list(
    objective = "binary:logistic",
    booster = "gbtree",
    eta = 0.01,
    max_depth = 10,
    min_child_weight = 20,
    subsample = 1,
    colsample_bytree = 1,
    lambda = 1,
    alpha = 0.5,
    eval_metric = "auc",
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
  
  l2g <- predict(xgb_model, dpred)
  
  ind=which(l2g>=0.5)
  ll=pred$studyLocusId[ind]
  ll1=table(ll)
  print(paste("N_CS with 1 gene:",sum(ll1==1)))
  print(paste("N_CS with 2 genes:",sum(ll1==2)))
  print(paste("N_CS with 2 genes:",sum(ll1>2)))
  
  probabilities <- predict(xgb_model, dtest)
  pr <- pr.curve(scores.class0 = probabilities, weights.class0 = y_test, curve = TRUE)
  
  auc_pr <- pr$auc.integral  # AUC-PR score
  print(paste("AUC-PR:", auc_pr))
  
  auc_result <- roc(y_test, probabilities)
  print(auc_result)
  
  matched_df <- merge(priv_l2g, cbind(pred[,c("studyLocusId","geneId")],l2g_new=l2g), by = c("studyLocusId","geneId"))
  print(paste("cor with 24.12:",cor(matched_df$l2g_new,matched_df$score)))
  cat("\n")
  cat("\n")
}

plot(matched_df$l2g_new, matched_df$score, 
     col = rgb(0, 0, 1, alpha = 0.005),  # Blue color with transparency
     pch = 16,  # Filled circle points
     xlab = "l2g_new", 
     ylab = "score", 
     main = "Scatter Plot of l2g_new vs. score")
