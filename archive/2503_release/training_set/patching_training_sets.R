#Making the patch of training sets

setwd("Desktop/training_set/")
library(data.table)

training_set_list=c("training_2503-testrun-1_all_string_extended_EGL_variants.tsv",
                    "training_2503-testrun-1_all_string_liberal_variants.tsv",
                    "training_2503-testrun-1_all_string_liberal_studies.tsv",
                    "training_2503-testrun-1_all_string_liberal_ALL_no_filter.tsv",
                    "training_2503-testrun-1_v5_1_no_exposion.tsv",
                    "training_2503-testrun-1_old_gs.tsv")

for (ts in training_set_list){
  fm=fread(ts,data.table=F)
  
  print(ts)
  
  n1=nrow(fm)
  tmp1=fm[fm$GSP==1,]
  
  colnms=c("geneId","efo_terms","variantId","vepMaximum","vepMean")
  colnms_to_round=c("eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
                    "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum")
  tmp1[colnms_to_round] <- lapply(tmp1[colnms_to_round], round, digits = 2)
  tmp1=tmp1[!(duplicated(tmp1[,c(colnms,colnms_to_round)])),]
  cs_to_keep=unique(tmp1$studyLocusId)
  fm=fm[fm$studyLocusId%in%cs_to_keep,]
  
  cs_list=unique(fm$studyLocusId)
  out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))

  cs_list_2=cs_list[out<3]
  fm=fm[fm$studyLocusId%in%cs_list_2,]
  cs_list=unique(fm$studyLocusId)
  n2=nrow(fm)
  print(paste("Reduced by:",round(((n1-n2)/n1)*100,2),"%"))
  write.table(fm,file=paste0("patched_",ts),sep="\t",quote=F,col.names = T,row.names = F)
}

#### making training set

set.seed(7777)
N_test=90
#### 1
ts="patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv"
test=fread(ts,data.table=F)
test_cs=unique(test$studyLocusId)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
ind=sample(1:length(test_geneid_efo),N_test)
test_geneid_efo=test_geneid_efo[ind]
cs_to_include=unique(test$studyLocusId[test$geneid_efo%in%test_geneid_efo])
test=test[test$studyLocusId%in%cs_to_include,]
write.table(test,file=paste0("test_",ts),sep="\t",quote=F,col.names = T,row.names = F)

#### 2
ts="patched_training_2503-testrun-1_old_gs.tsv"
test=fread(ts,data.table=F)
test_cs=unique(test$studyLocusId)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
ind=sample(1:length(test_geneid_efo),N_test)
test_geneid_efo=test_geneid_efo[ind]
cs_to_include=unique(test$studyLocusId[test$geneid_efo%in%test_geneid_efo])
test=test[test$studyLocusId%in%cs_to_include,]
write.table(test,file=paste0("test_",ts),sep="\t",quote=F,col.names = T,row.names = F)

#### 3
ts="patched_training_2503-testrun-1_v5_1_no_exposion.tsv"
test=fread(ts,data.table=F)
test_cs=unique(test$studyLocusId)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])
ind=sample(1:length(test_geneid_efo),N_test)
test_geneid_efo=test_geneid_efo[ind]
cs_to_include=unique(test$studyLocusId[test$geneid_efo%in%test_geneid_efo])
test=test[test$studyLocusId%in%cs_to_include,]
write.table(test,file=paste0("test_",ts),sep="\t",quote=F,col.names = T,row.names = F)

##### DESCRIPTIVES

training_set_list=c("patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_variants.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_studies.tsv",
                    "patched_training_2503-testrun-1_all_string_liberal_ALL_no_filter.tsv",
                    "patched_training_2503-testrun-1_v5_1_no_exposion.tsv",
                    "patched_training_2503-testrun-1_old_gs.tsv")
test_set="test_patched_training_2503-testrun-1_all_string_extended_EGL_variants.tsv"

#### descriptives

for (ts in c(training_set_list,test_set)){
  fm=fread(ts,data.table=F)
  
  print(ts)
  print(paste("GSN:",table(fm$GSP)[1]))
  print(paste("GSP:",table(fm$GSP)[2]))
  print(paste("class imbalance: 1 :",table(fm$GSP)[1]/table(fm$GSP)[2]))
  print(paste("unique geneId_efo pairs for GSP:",length(unique((fm$geneid_efo[fm$GSP==1])))))
  print(paste("unique studyLocusId:",length(unique((fm$studyLocusId)))))
  
  out=NULL
  gidefo=unique((fm$geneid_efo[fm$GSP==1]))
  for (gid in gidefo){
    x <- fm[fm$geneid_efo == gid & fm$GSP==1, ]
    out=c(out,sum(x$distanceSentinelFootprintNeighbourhood==1 | x$distanceSentinelTssNeighbourhood==1)>0)
  }
  
  ll=table(out)
  print(paste("The proportion of nearest to footrpint or tss gene for GSP:",ll[2]/(ll[1]+ll[2])))
  cat("\n")
  cat("\n")
}

#### descriptives excluding test set
test=fread(test_set,data.table=F)
test_cs=unique(test$studyLocusId)
test_geneid_efo=unique(test$geneid_efo[test$GSP==1])

for (ts in c(training_set_list)){
  fm=fread(ts,data.table=F)
  
  fm1=fm[fm$GSP==1,]
  fm1=fm1[!(fm1$geneid_efo%in%test_geneid_efo),]
  fm1_cs=unique(fm1$studyLocusId)
  fm=fm[fm$studyLocusId%in%fm1_cs,]
  
  print(ts)
  print(paste("GSN:",table(fm$GSP)[1]))
  print(paste("GSP:",table(fm$GSP)[2]))
  print(paste("class imbalance: 1 :",table(fm$GSP)[1]/table(fm$GSP)[2]))
  print(paste("unique geneId_efo pairs for GSP:",length(unique((fm$geneid_efo[fm$GSP==1])))))
  print(paste("unique studyLocusId:",length(unique((fm$studyLocusId)))))
  
  out=NULL
  gidefo=unique((fm$geneid_efo[fm$GSP==1]))
  for (gid in gidefo){
    x <- fm[fm$geneid_efo == gid & fm$GSP==1, ]
    out=c(out,sum(x$distanceSentinelFootprintNeighbourhood==1 | x$distanceSentinelTssNeighbourhood==1)>0)
  }
  
  ll=table(out)
  print(paste("The proportion of nearest to footrpint or tss gene for GSP:",ll[2]/(ll[1]+ll[2])))
  cat("\n")
  cat("\n")
}





