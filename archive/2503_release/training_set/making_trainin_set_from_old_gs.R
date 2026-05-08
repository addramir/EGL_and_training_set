setwd("Desktop/training_set/")
library(data.table)

x=fread("training_2503-testrun-1_all_string_liberal_ALL_no_filter.tsv",data.table = F)

x=cbind(x,id1=paste(x$geneId,x$studyId,x$variantId,sep="_"))
table(duplicated(x$id1))

old_gs=fread("old_gs_high_medium_GSP_only.csv/part-00000-6ec76868-71cd-4e85-b8e7-9eb595be631a-c000.csv",data.table=FALSE)
old_gs=cbind(old_gs,variantId=paste(old_gs$chromosome,old_gs$position,old_gs$reference,old_gs$alternative,sep="_"))

table(old_gs$variantId%in%x$variantId)

old_gs=cbind(old_gs,id1=paste(old_gs$geneId,old_gs$otg_id,old_gs$variantId,sep="_"))
old_gs=cbind(old_gs,id2=paste(old_gs$geneId,old_gs$gwas_catalog_id,old_gs$variantId,sep="_"))

table(duplicated(old_gs$id2))
old_gs=old_gs[-which(duplicated(old_gs$id2)),]
table(old_gs$id2%in%x$id1)

old_gs=old_gs[old_gs$id2%in%x$id1,]
ind=match(old_gs$id2,x$id1)
table(old_gs$id2==x$id1[ind])
old_gs=cbind(old_gs,x[ind,-which(colnames(x)%in%colnames(old_gs))])

table(old_gs$GSP)
old_gs=old_gs[old_gs$GSP==1,]

old_gs_gsp_only=old_gs
length(unique((old_gs$geneid_efo)))
length(unique((old_gs$studyLocusId)))

cs_to_include=unique(old_gs$studyLocusId)
fm=x[x$studyLocusId%in%cs_to_include,]
ind=which(!(fm$id1%in%old_gs$id2) & fm$GSP==1)
fm=fm[-ind,]

table(fm$GSP)
table(old_gs$GSP)

cs_list=unique(fm$studyLocusId)
out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
table(out)

cs_list_2=cs_list[out<3]
fm=fm[fm$studyLocusId%in%cs_list_2,]
cs_list=unique(fm$studyLocusId)

boxplot(fm$sum_col~fm$GSP)

out <- t(sapply(cs_list, function(cs) {
  x <- fm[fm$studyLocusId == cs, ]
  
  c(
    GSP = x$geneId[which(x$GSP == 1)[1]],
    nearest_tss = x$geneId[which.max(x$distanceSentinelTss)[1]],
    nearest_footprint = x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  )
}))

table(out[,1]==out[,2] | out[,1]==out[,3])

write.table(fm,file="training_2503-testrun-1_old_gs.tsv",sep="\t",quote=F,col.names = T,row.names = F)



