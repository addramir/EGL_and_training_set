setwd("Desktop/training_set/")

library(data.table)

#LOAD AND PREPROCESSING
gsp=fread("2503-testrun-1-liberal_GSP_only.csv/part-00000-dafe1aac-b3b9-4a65-b778-ed5298c59c3f-c000.csv",
          data.table=F, header = T,sep="\t",quote="")
colnames(gsp)[which(colnames(gsp)=="diseaseId")]="efo_terms"

gsp=gsp[gsp[,"isProteinCoding"]==1,]

cs_list=unique(gsp$studyLocusId)

fm=fread("2503-testrun-1-liberal_FM.csv/part-00000-7a454abe-88d7-4a7f-9afb-5134bf889d71-c000.csv",
         data.table=F, header = T,sep="\t",quote="")
fm=fm[fm[,"isProteinCoding"]==1,]

fm=fm[fm$studyLocusId%in%cs_list,]

gsp_efo=gsp[!duplicated(gsp$studyLocusId),][,c("studyLocusId","efo_terms","beta","pValueExponent")]
ind=match(fm$studyLocusId,gsp_efo$studyLocusId)
table(gsp_efo$studyLocusId[ind]==fm$studyLocusId)
fm=cbind(fm,efo_terms=gsp_efo$efo_terms[ind])
#fm=cbind(fm,beta=abs(gsp_efo$beta[ind]))
#fm=cbind(fm,pval=abs(gsp_efo$pValueExponent[ind]))
fm=cbind(fm,"geneid_efo"=paste0(fm$geneId,"_",fm$efo_terms))

GSP=paste0(gsp$geneId,"_",gsp$efo_terms)
GSP=unique(GSP)

fm=cbind(fm,GSP=0)
fm[fm$geneid_efo%in%GSP,"GSP"]=1
cs_list=unique(fm$studyLocusId)

#ADDING studyId
#dim(gsp)
#gsp1=gsp[!duplicated(gsp$studyLocusId),]
#dim(gsp1)
#ind=match(fm$studyLocusId,gsp1$studyLocusId)
#table(fm$studyLocusId==gsp1$studyLocusId[ind])
#fm$studyId=gsp1$studyId[ind]

# FILTERING BY SUMCOL
fm[fm[,"vepMaximum"]<=0.1,"vepMaximum"]=0
x=c("distanceSentinelFootprint","eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
    "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum","vepMaximum")
ll=apply(fm[,x],MAR=1,sum)
fm=cbind(fm,sum_col=ll)

fm=fm[-which(fm$GSP==1 & fm$sum_col<0.95),]
cs_list=unique(fm$studyLocusId)
out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
table(out)
#0     1     2     3 
#2147 12161   380    95 

cs_list_2=cs_list[out>0]
fm=fm[fm$studyLocusId%in%cs_list_2,]
cs_list=unique(fm$studyLocusId)

boxplot(fm$sum_col~fm$GSP)

### CHECKING SINGLE STUDIES for GENE_EFO
fm1=fm[fm$GSP==1,]
geneid_efo_list=unique(fm1$geneid_efo)
ll <- sapply(geneid_efo_list, function(geidefo) length(unique(fm1$studyId[fm1$geneid_efo == geidefo])))
length(geneid_efo_list)
table(ll==1)


single_geneid_efo=names(ll)[ll==1]

fm2=fm1[fm1$geneid_efo%in%single_geneid_efo,]
cs_signle_geneid_efo=fm2$studyLocusId
cs_list_3=unique(cs_signle_geneid_efo)
fm2=fm[fm$studyLocusId%in%cs_list_3,]

#fm3=fm2[fm2$GSP==1,]
#cs_list_3=unique(fm3$studyLocusId)
#fm2=fm[fm$studyLocusId%in%cs_list_3,]

out <- t(sapply(cs_list_3, function(cs) {
  x <- fm2[fm2$studyLocusId == cs, ]
  
  c(
    GSP = x$geneId[which(x$GSP == 1)[1]],
    nearest_tss = x$geneId[which.max(x$distanceSentinelTss)[1]],
    nearest_footprint = x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  )
}))

table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#240   131


not_single_geneid_efo=names(ll)[ll>1]
fm1=fm[fm$GSP==1,]
fm2=fm1[fm1$geneid_efo%in%not_single_geneid_efo,]
cs_not_signle_geneid_efo=fm2$studyLocusId
cs_list_3=unique(cs_not_signle_geneid_efo)
fm2=fm[fm$studyLocusId%in%cs_list_3,]

fm4=fm2[fm2$GSP==1,]
geneid_efo_list=unique(fm4$geneid_efo)
ll <- sapply(geneid_efo_list, function(geidefo) length(unique(fm4$studyId[fm4$geneid_efo == geidefo])))
length(geneid_efo_list)
table(ll==1)



#fm3=fm2[fm2$GSP==1,]
#cs_list_3=unique(fm3$studyLocusId)
#fm2=fm[fm$studyLocusId%in%cs_list_3,]

out <- t(sapply(cs_list_3, function(cs) {
  x <- fm2[fm2$studyLocusId == cs, ]
  
  c(
    GSP = x$geneId[which(x$GSP == 1)[1]],
    nearest_tss = x$geneId[which.max(x$distanceSentinelTss)[1]],
    nearest_footprint = x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  )
}))

table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#3871  8398 


### FILTER by string
fm3=fm2
cs_list_4=unique(fm3$studyLocusId)
cs=cs_list_4[1]
out <- vapply(cs_list_4, function(cs) sum(fm3$GSP[fm3$studyLocusId == cs] == 1), numeric(1))

fm4=fm3[fm3$GSP==1,]
geneid_efo_list=unique(fm4$geneid_efo)
ll <- sapply(geneid_efo_list, function(geidefo) length(unique(fm4$studyId[fm4$geneid_efo == geidefo])))
length(geneid_efo_list)
table(ll==1)

table(out)
#    1     2     3 
#11814   364    91 


strng=fread("interactions.tsv/part-00000-74b0764c-70d6-4851-aba5-01c2cc225866-c000.csv",data.table=F)

dim(fm3)
#[1] 193559     39  
table(fm3$GSP)
#0      1 
#180744  12815  
i=1
for (cs in cs_list_4){
  x=fm3[fm3$studyLocusId==cs,]
  y=x[x$GSP==1,"geneId"]
  
  to_exclude=NULL
  ind=which(strng$targetA%in%y)
  if (length(ind)>0){
    to_exclude=unique(strng$targetB[ind])
  }
  
  ind=which(strng$targetB%in%y)
  if (length(ind)>0){
    to_exclude=unique(c(to_exclude,strng$targetA[ind]))
  }
  
  if(length(to_exclude)>0){
    ind=which(to_exclude%in%y)
    if (length(ind)>0){
      to_exclude=to_exclude[-ind]
    }
    ind=which(fm3$studyLocusId==cs & fm3$geneId%in%to_exclude)
    if (length(ind)>0){
      fm3=fm3[-ind,]
    }
  }
}
dim(fm3)
#[1] 189197     39
length(unique(fm3$geneid_efo[fm3$GSP==1]))
#494

write.table(fm3,file="training_2503-testrun-1_all_string_liberal_studies.tsv",sep="\t",quote=F,col.names = T,row.names = F)


