setwd("~/Projects/gentropy/notebooks/2506_release/gwas_catalog_curation_v2/")

library(data.table)

x=fread("20250612_output_from_curation.txt",data.table=F)

#example=fread("GWAS_Catalog_study_curation.tsv",data.table=F)

clnms=c("studyId", "studyType","analysisFlag","qualityControl","isCurated",
        "pubmedId","publicationTitle","traitFromSource")

colnames(x)

colnames(x)[colnames(x)=="PUBMED ID"]="pubmedId"
colnames(x)[colnames(x)=="STUDY"]="publicationTitle"
colnames(x)[colnames(x)=="STUDY ACCESSION"]="studyId"
colnames(x)[colnames(x)=="DISEASE/TRAIT"]="traitFromSource"

x=cbind(x,qualityControl=NA)

x=x[,clnms]

fwrite(x,"20250612_output_curation_formated.tsv",sep="\t",quote=F)
