setwd("Desktop/training_set/")

library(data.table)
library(rjson)

full_fm=fread("training_2503-testrun-1_all_string_liberal_egl.tsv",data.table=FALSE)
# OLD MODEL
fm=fromJSON(file="WandB/old_model/conservative_max200/featureMatrix_13_42ecdcb6bbe329f8158d.table.json",simplify = TRUE)
l=unlist(fm$data)
l1=t(matrix(l,nrow=length(fm$columns)))
ts=as.data.frame(l1)
colnames(ts)=fm$columns
ts=cbind(ts,geneid_efo=paste(ts$geneId,ts$traitFromSourceMappedId,sep="_"))
length(unique(ts$geneid_efo[ts$goldStandardSet=="positive"]))
length(unique(ts$studyLocusId))

table(unique(ts$studyLocusId)%in%unique(full_fm$studyLocusId))
