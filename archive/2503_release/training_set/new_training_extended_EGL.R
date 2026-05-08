setwd("Desktop/training_set/")

library(data.table)

#LOAD AND PREPROCESSING
gsp=fread("2503-testrun-1-extended_EGL_GSP_only.csv/part-00000-b55ffdee-9336-42c6-b7b2-f7540f537cec-c000.csv",
          data.table=F, header = T,sep="\t",quote="")
colnames(gsp)[which(colnames(gsp)=="diseaseId")]="efo_terms"

gsp=gsp[gsp[,"isProteinCoding"]==1,]

cs_list=unique(gsp$studyLocusId)

fm=fread("2503-testrun-1-extended_EGL_FM.csv/part-00000-674d3cd4-8972-4d92-b1dd-d00b7a010c67-c000.csv",
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

#filtering by dist 0.1
fm=fm[-which(fm$GSP==1 & fm$distanceSentinelFootprint<0.1),]
cs_list=unique(fm$studyLocusId)
out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
table(out)
#0     1     2     3     4 
#298 13850   503   151     7 

cs_list_2=cs_list[out>0]
fm=fm[fm$studyLocusId%in%cs_list_2,]
cs_list=unique(fm$studyLocusId)


# dupl by variat
fm1=fm[fm$GSP==1,]
ll=paste(fm1$variantId,fm1$geneId,fm1$efo_terms,sep="_")
dupl_ll=unique(ll[which(duplicated(ll))])
fm1=fm1[ll%in%dupl_ll,]
length(unique(fm1$geneid_efo))
# 442

cs_to_include=unique(fm1$studyLocusId)
fm=fm[fm$studyLocusId%in%cs_to_include,]
cs_list=unique(fm$studyLocusId)
out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
table(out)
#1     2     3     4     5 
#10586   305    91     3    17 

cs_list_2=cs_list[out<3]
fm=fm[fm$studyLocusId%in%cs_list_2,]
cs_list=unique(fm$studyLocusId)
out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
table(out)

#sumcol
fm[fm[,"vepMaximum"]<=0.1,"vepMaximum"]=0
x=c("distanceSentinelFootprint","eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
    "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum","vepMaximum")
ll=apply(fm[,x],MAR=1,sum)
fm=cbind(fm,sum_col=ll)

#fm=fm[-which(fm$GSP==1 & fm$sum_col<0.95),]
#cs_list=unique(fm$studyLocusId)
#out <- vapply(cs_list, function(cs) sum(fm$GSP[fm$studyLocusId == cs] == 1), numeric(1))
#table(out)
#0     1     2     3 
#2147 12161   380    95 

#cs_list_2=cs_list[out>0]
#fm=fm[fm$studyLocusId%in%cs_list_2,]
#cs_list=unique(fm$studyLocusId)

boxplot(fm$sum_col~fm$GSP)

### CHECKING SINGLE STUDIES for GENE_EFO

out <- t(sapply(cs_list, function(cs) {
  x <- fm[fm$studyLocusId == cs, ]
  
  c(
    GSP = x$geneId[which(x$GSP == 1)[1]],
    nearest_tss = x$geneId[which.max(x$distanceSentinelTss)[1]],
    nearest_footprint = x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  )
}))

table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#3624  7267 


### FILTER by string
fm3=fm
cs_list_4=unique(fm3$studyLocusId)
cs=cs_list_4[1]
out <- vapply(cs_list_4, function(cs) sum(fm3$GSP[fm3$studyLocusId == cs] == 1), numeric(1))

fm4=fm3[fm3$GSP==1,]
geneid_efo_list=unique(fm4$geneid_efo)
ll <- sapply(geneid_efo_list, function(geidefo) length(unique(fm4$studyId[fm4$geneid_efo == geidefo])))
length(geneid_efo_list)
#485
table(ll==1)

table(out)
#1     2 
#10586   305 


#strng=fread("interactions.tsv/part-00000-74b0764c-70d6-4851-aba5-01c2cc225866-c000.csv",data.table=F)
strng=fread("interactions_005.csv/part-00000-b3e13269-bb08-4c9b-8469-8f7b9a4c9e6c-c000.csv",data.table=F)

dim(fm3)
#173267     39
table(fm3$GSP)
#     0      1 
#162071  11196  
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
#169085     39

#128658     39 - new on 0.05
length(unique(fm3$geneid_efo[fm3$GSP==1]))
#485

write.table(fm3,file="training_2503-testrun-1_all_string_005_extended_EGL_variants.tsv",sep="\t",quote=F,col.names = T,row.names = F)




