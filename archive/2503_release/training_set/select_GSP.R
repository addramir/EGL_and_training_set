setwd("Desktop/training_set/")

library(data.table)
gsp=fread("2503-testrun-1_GSP_only.csv/part-00000-552d8dbc-d875-436c-9897-c09becab1e3e-c000.csv",
          data.table=F, header = T,sep="\t",quote="")
colnames(gsp)[which(colnames(gsp)=="diseaseId")]="efo_terms"

gsp=gsp[gsp[,"isProteinCoding"]==1,]

cs_list=unique(gsp$studyLocusId)

fm=fread("2503-testrun-1_FM.csv/part-00000-ef160c55-a059-4972-a4fa-ac1263e0c5aa-c000.csv",
         data.table=F, header = T,sep="\t",quote="")
fm=fm[fm[,"isProteinCoding"]==1,]

fm=fm[fm$studyLocusId%in%cs_list,]

gsp_efo=gsp[!duplicated(gsp$studyLocusId),][,c("studyLocusId","efo_terms","beta","pValueExponent")]
ind=match(fm$studyLocusId,gsp_efo$studyLocusId)
table(gsp_efo$studyLocusId[ind]==fm$studyLocusId)
fm=cbind(fm,efo_terms=gsp_efo$efo_terms[ind])
fm=cbind(fm,beta=abs(gsp_efo$beta[ind]))
fm=cbind(fm,pval=abs(gsp_efo$pValueExponent[ind]))
fm=cbind(fm,"geneid_efo"=paste0(fm$geneId,"_",fm$efo_terms))

GSP=paste0(gsp$geneId,"_",gsp$efo_terms)
GSP=unique(GSP)

fm=cbind(fm,GSP=0)
fm[fm$geneid_efo%in%GSP,"GSP"]=1
cs_list=unique(fm$studyLocusId)

fm=fm[-which(fm$GSP==1 & fm$distanceSentinelFootprint<0.5),]
cs_list=unique(fm$studyLocusId)

cs=cs_list[1]
out=NULL
for (cs in cs_list){
  x=fm[fm$studyLocusId==cs,]
  out=c(out,sum(x$GSP==1))
}
table(out)
#0     1     2     3     4 
#301 13847   503   151     7 


cs_list_2=cs_list[out>0]
fm=fm[fm$studyLocusId%in%cs_list_2,]
cs_list=unique(fm$studyLocusId)


out=array(NA,c(length(cs_list),3))
colnames(out)=c("GSP","nearest_tss","nearest_footprint")
i=1
for (cs in cs_list){
  x=fm[fm$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  i=i+1
}

table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#5985  8523 


#write.table(fm,file="training_v10_6_full_full.tsv",sep="\t",quote=F,col.names = T,row.names = F)
#fm2=fread("training_v10_3_full.tsv",sep="\t",data.table = F)


fm1=fm[fm$GSP==1,]
length(unique(fm1$geneid_efo))
l=table(fm1$geneid_efo)
hist(l)
table(l==1)
#FALSE  TRUE 
#673   473 

single_geneid_efo=names(l)[l==1]
not_single_geneid_efo=names(l)[l>1]

fm2=fm1[fm1$geneid_efo%in%single_geneid_efo,]
cs_signle_geneid_efo=fm2$studyLocusId
cs_list_3=unique(cs_signle_geneid_efo)
fm2=fm[fm$studyLocusId%in%cs_list_3,]

fm3=fm2[fm2$GSP==1,]
#fm3=fm3[fm3$sum_col>=1.5,]
cs_list_3=unique(fm3$studyLocusId)
fm2=fm[fm$studyLocusId%in%cs_list_3,]

out=array(NA,c(length(cs_list_3),3))
colnames(out)=c("GSP","nearest_tss","nearest_footprint")
i=1
for (cs in cs_list_3){
  x=fm2[fm2$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  i=i+1
}
table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#338    93 

not_single_geneid_efo=names(l)[l>1]
fm2=fm1[fm1$geneid_efo%in%not_single_geneid_efo,]
cs_not_signle_geneid_efo=fm2$studyLocusId
cs_list_3=unique(cs_not_signle_geneid_efo)
fm2=fm[fm$studyLocusId%in%cs_list_3,]
cs_list_3=unique(fm2$studyLocusId)
out=array(NA,c(length(cs_list_3),3))
colnames(out)=c("GSP","nearest_tss","nearest_footprint")
i=1
for (cs in cs_list_3){
  x=fm2[fm2$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  i=i+1
}
table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#5651  8432 

#write.table(fm2,file="training_v24106_full_no_single_cs.tsv",sep="\t",quote=F,col.names = T,row.names = F)
#fm2=fread("training_v10_6_full_no_single_cs.tsv",sep="\t",data.table = F)

#######
fm2[fm2[,"vepMaximum"]<=0.1,"vepMaximum"]=0
x=c("distanceSentinelFootprint","eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
   "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum","vepMaximum")
ll=apply(fm2[,x],MAR=1,sum)
fm2=cbind(fm2,sum_col=ll)
x1=fm2$sum_col[fm2$GSP==1]
x2=fm2$sum_col[fm2$GSP==0]
summary(x2)
ll=cbind(c(x1,x2),c(rep(1,length(x1)),rep(0,length(x2))))
boxplot(ll)

out=array(NA,c(length(cs_list_3),4))
colnames(out)=c("GSP","nearest_tss","nearest_footprint","sum_col")
i=1
for (cs in cs_list_3){
  x=fm2[fm2$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  out[i,4]=x$geneId[which.max(x$sum_col)[1]]
  i=i+1
}

table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#5651  8432 

table(out[,1]==out[,4] | out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#5043  9040 

#select CS _ sum_col

fm2=fm2[-which(fm2$GSP==1 & fm2$sum_col<0.95),]
cs_list=unique(fm2$studyLocusId)

cs=cs_list[1]
out=NULL
for (cs in cs_list){
  x=fm2[fm2$studyLocusId==cs,]
  out=c(out,sum(x$GSP==1))
}
table(out)
#0     1     2     3 
#1681 11919   367    91


cs_list_2=cs_list[out>0]
fm2=fm2[fm2$studyLocusId%in%cs_list_2,]
cs_list=unique(fm2$studyLocusId)


### FILTER by string
fm3=fm2
cs_list_4=unique(fm3$studyLocusId)
cs=cs_list_4[1]
out=NULL
for (cs in cs_list_4){
  x=fm3[fm3$studyLocusId==cs,]
  out=c(out,sum(x$GSP==1))
}
table(out)
#1     2     3 
#11919   367    91 


strng=fread("interactions.tsv/part-00000-74b0764c-70d6-4851-aba5-01c2cc225866-c000.csv",data.table=F)

dim(fm3)
#[1] 195440     36
table(fm3$GSP)
#0      1 
#182514  12926  
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
#[1] 191062     36
length(unique(fm3$geneid_efo[fm3$GSP==1]))
#572

write.table(fm3,file="training_2503-testrun-1_all_string_conservative_egl.tsv",sep="\t",quote=F,col.names = T,row.names = F)

# SELECT first 200
nmax=200
x=fread("training_2503-testrun-1_all_string_liberal_egl.tsv",data.table=FALSE)
y=x[x$GSP==1,]
l=table(y$geneid_efo)
geidlist=names(l[l>nmax])

geid=geidlist[1]
cs_to_exclude=NULL
for (geid in geidlist){
  y1=y[y$geneid_efo==geid,]
  y1=y1[order(y1$sum_col,decreasing = TRUE),]
  y2=y1[-(1:nmax),]
  cs_to_exclude=c(cs_to_exclude,y2$studyLocusId)
}

cs_to_exclude=unique(cs_to_exclude)

x1=x[!x$studyLocusId%in%cs_to_exclude,]
write.table(x1,file="training_2503-testrun-1_all_string_liberal_egl_max200.tsv",sep="\t",quote=F,col.names = T,row.names = F)


#STOPPPPPPP


#IGNORE
fm3=fm2[fm2$GSP==1,]
egl=unique(fm3$geneid_efo)
cs_to_use=NULL
for (gene in egl){
  x=fm3[fm3$geneid_efo==gene,]
  l=x$studyLocusId[which.max(x$sum_col)[1]]
  cs_to_use=c(cs_to_use,l)
}

fm3=fm2[fm2$studyLocusId%in%cs_to_use,]
cs_list_4=unique(fm3$studyLocusId)
out=array(NA,c(length(cs_list_4),4))
colnames(out)=c("GSP","nearest_tss","nearest_footprint","sum_col")
i=1
for (cs in cs_list_4){
  x=fm3[fm3$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  out[i,4]=x$geneId[which.max(x$sum_col)[1]]
  i=i+1
}
table(out[,1]==out[,2])
table(out[,1]==out[,3])
table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#173   259 
table(out[,1]==out[,4])
#FALSE  TRUE 
#189   243
table(out[,1]==out[,4] | out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#157   275 

write.table(fm2,file="training_v10_v6_sum_col_best_cs_selected.tsv",sep="\t",quote=F,col.names = T,row.names = F)

fm3=fread("training_v10_v6_sum_col_best_cs_selected.tsv",data.table=F)
#STOP IGNORE



#select CS - dist
fm3=fm2[fm2$GSP==1,]
egl=unique(fm3$geneid_efo)
cs_to_use=NULL
for (gene in egl){
  x=fm3[fm3$geneid_efo==gene,]
  #l=x$studyLocusId[which.max(x$pval)[1]]
  y=x[order(x$distanceSentinelFootprint,x$sum_col,decreasing = T),]
  l=y$studyLocusId[1]
  cs_to_use=c(cs_to_use,l)
}

fm3=fm2[fm2$studyLocusId%in%cs_to_use,]
cs_list_4=unique(fm3$studyLocusId)
out=array(NA,c(length(cs_list_4),4))
colnames(out)=c("GSP","nearest_tss","nearest_footprint","sum_col")
i=1
for (cs in cs_list_4){
  x=fm3[fm3$studyLocusId==cs,]
  out[i,1]=x$geneId[which(x$GSP==1)[1]]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  out[i,4]=x$geneId[which.max(x$sum_col)[1]]
  i=i+1
}
table(out[,1]==out[,2])
table(out[,1]==out[,3])
table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#271   290
table(out[,1]==out[,4])
#FALSE  TRUE 
#305   256 
table(out[,1]==out[,4] | out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#263   298 
write.table(fm3,file="training_v2_dist_sum.tsv",sep="\t",quote=F,col.names = T,row.names = F)




#select CS - dist - pval
fm3=fm2[fm2$GSP==1,]
egl=unique(fm3$geneid_efo)
cs_to_use=NULL
for (gene in egl){
  x=fm3[fm3$geneid_efo==gene,]
  #l=x$studyLocusId[which.max(x$pval)[1]]
  y=x[order(x$distanceSentinelFootprint,x$pval,x$sum_col,decreasing = T),]
  l=y$studyLocusId[1]
  cs_to_use=c(cs_to_use,l)
}

fm3=fm2[fm2$studyLocusId%in%cs_to_use,]
cs_list_4=unique(fm3$studyLocusId)
out=array(NA,c(length(cs_list_4),4))
colnames(out)=c("GSP","nearest_tss","nearest_footprint","sum_col")
i=1
for (cs in cs_list_4){
  x=fm3[fm3$studyLocusId==cs,]
  out[i,1]=x$geneId[x$GSP==1]
  out[i,2]=x$geneId[which.max(x$distanceSentinelTss)[1]]
  out[i,3]=x$geneId[which.max(x$distanceSentinelFootprint)[1]]
  out[i,4]=x$geneId[which.max(x$sum_col)[1]]
  i=i+1
}
table(out[,1]==out[,2])
table(out[,1]==out[,3])
table(out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#202   250 
table(out[,1]==out[,4])
#FALSE  TRUE 
#215   221 
table(out[,1]==out[,4] | out[,1]==out[,2] | out[,1]==out[,3])
#FALSE  TRUE 
#198   253 
write.table(fm3,file="training_v3_dist_pval_sum.tsv",sep="\t",quote=F,col.names = T,row.names = F)


x=c("distanceSentinelFootprint","eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
    "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum","vepMaximum")
l=princomp(fm3[,x])
plot(l$scores)
points(l$scores[fm3$GSP==1,],col="red")




x1=fread("../training_v1_sum_col.tsv",data.table=F)
x2=fread("../training_v2_dist_sum.tsv",data.table=F)
x3=fread("../training_v2_dist_pval_sum.tsv",data.table=F)

x1c=unique(x1$studyLocusId)
x2c=unique(x2$studyLocusId)
x3c=unique(x3$studyLocusId)
library(VennDiagram)
venn.diagram(x = list(x1c,x2c,x3c),filename="tmp.png",category.names = c("x 1" , "x 2 " , "x 3"))


install.packages("jsonlite")
library(jsonlite)

# Read the JSON file (using the correct file path)
json_data <- fromJSON("../featureMatrix_12_72e1caafae1514639eeb.table.json",flatten=TRUE)
df=json_data[[2]]
colnames(df)=json_data[[1]]
df=as.data.frame(df)
dim(df)
table(df$goldStandardSet)
length(unique(df$studyLocusId))
cs_list_old=unique(df$studyLocusId)
cs=cs_list_old[1]
out=NULL
for (cs in cs_list_old){
  x=df[df$studyLocusId==cs,]
  out=c(out,sum(x$goldStandardSet=="positive"))
}
table(out)
cs_list_old=cs_list_old[out==1]
df=df[df$studyLocusId%in%cs_list_old,]

fm2=df
fm2[fm2[,"vepMaximum"]<=0.1,"vepMaximum"]=0
x=c("distanceSentinelFootprint","eQtlColocClppMaximum","pQtlColocClppMaximum","sQtlColocClppMaximum",
    "eQtlColocH4Maximum","pQtlColocH4Maximum","sQtlColocH4Maximum","vepMaximum")

l1=fm2[,x]
l1=apply(l1,MAR=2,as.numeric)
ll=apply(l1,MAR=1,sum)
fm2=cbind(fm2,sum_col=ll)
x1=fm2$sum_col[fm2$goldStandardSet=="positive"]
x2=fm2$sum_col[fm2$goldStandardSet=="negative"]
summary(x2)
ll=cbind(c(x1,x2),c(rep(1,length(x1)),rep(0,length(x2))))
boxplot(ll)



