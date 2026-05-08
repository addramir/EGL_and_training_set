gsutil ls 'gs://gwas_catalog_inputs/raw_summary_statistics/**/*h.tsv.gz' > gsutil_list.txt
gsutil cp gs://gwas_catalog_inputs/gwas_catalog_download_studies.tsv .
gsutil cp gs://genetics-portal-dev-analysis/yt4/20241004_gwascat_curation.tsv .

#gsutil cp 20241004_input_for_curation.tsv gs://genetics-portal-dev-analysis/yt4/20241004_input_for_curation.tsv

R
library(data.table)
setwd("Desktop/gwas_cat_harmon/")
gwascat=fread("gwas-catalog-v1.0.3.1-studies-r2024-12-19.tsv",data.table=FALSE,quote="")

lss=fread("gsutil_list.txt",data.table=FALSE,header=FALSE)
lss=lss[,1]

l1=gsub(pattern="gs://gwas_catalog_inputs/raw_summary_statistics/",repl="",x=lss)

ll1=strsplit(l1,split="/")

l2=c()
for (i in 1:length(ll1)){
	x=ll1[[i]]
	l2=c(l2,x[2])
}

table(l2%in%gwascat[,"STUDY ACCESSION"])
#FALSE  TRUE 
#8267 72643

for_curation=gwascat[gwascat[,"STUDY ACCESSION"]%in%l2,]

old_cur=fread("yt4_gwas_catalog_curation_20241120_output_curation.tsv",data.table=FALSE,quote="")

old_cur=old_cur[old_cur[,"isCurated"]==TRUE,]

old_cur=old_cur[old_cur$studyId%in%for_curation[,"STUDY ACCESSION"],]

ind=match(old_cur$studyId,for_curation[,"STUDY ACCESSION"])
table(old_cur$studyId==for_curation[ind,"STUDY ACCESSION"])

#studyType analysisFlag qualityControl isCurated

for_curation[,"studyType"]=""
for_curation[,"analysisFlag"]=""
for_curation[,"isCurated"]=FALSE

for_curation[ind,"studyType"]=old_cur[,"studyType"]
for_curation[ind,"analysisFlag"]=old_cur[,"analysisFlag"]
for_curation[ind,"isCurated"]=TRUE

clnms=c("PUBMED ID","STUDY","DISEASE/TRAIT","STUDY ACCESSION","FULL SUMMARY STATISTICS","studyType","analysisFlag","isCurated")


#x=fread("CUR.tsv",data.table=F)
x=for_curation[,clnms]

unique_pmid=unique(x$`PUBMED ID`)

for (upim in unique_pmid){
  ind=which(x$`PUBMED ID`==upim)
  
  if (sum(x[ind,"isCurated"])>0){
    x[ind,"isCurated"]=TRUE
  }
  
  ll=unique(x[ind,"studyType"])
  ll=ll[ll!=""]
  if (length(ll)>0){
    x[ind,"studyType"]=ll[1]
    x[ind,"isCurated"]=TRUE
  }
  
  ll=unique(x[ind,"analysisFlag"])
  ll=ll[ll!=""]
  if (length(ll)>0){
    x[ind,"analysisFlag"]=ll[1]
    x[ind,"isCurated"]=TRUE
  }
  
  
}

dim(for_curation)
#[1] 72643    28

table(for_curation$isCurated)
#FALSE  TRUE 
#14792 57851 


fwrite(x=x,file="20241219_input_for_curation.tsv",sep="\t",quote=FALSE)
