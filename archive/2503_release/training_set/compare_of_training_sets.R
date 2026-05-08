library(data.table)
x1=fread("training_v24106_sum_col_all_cs_string.tsv",data.table=F)
x2=fread("training_v2502_sum_col_all_cs_string.tsv",data.table=F)

l1=paste(x1$studyLocusId,x1$geneId,sep="_")
l2=paste(x2$studyLocusId,x2$geneId,sep="_")
table(l1%in%l2)

table(x1$GSP)
table(x2$GSP)

cs1=unique(x1$studyLocusId)
cs2=unique(x2$studyLocusId)

table(cs1%in%cs2)

x11=x1[x1$GSP==1,]
x21=x2[x2$GSP==1,]

gcp1=unique(paste(x11$geneId,x11$efo_terms,sep="_"))
gcp2=unique(paste(x21$geneId,x21$efo_terms,sep="_"))

table(gcp1%in%gcp2)

table(gcp2%in%gcp1)
