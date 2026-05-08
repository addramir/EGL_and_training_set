setwd("~/Projects/gentropy/notebooks/gentropy_paper/")

library(data.table)
library(arrow)
#x=fread("data/protein_coding_eqtls_l2g_fm_based_rare_varint_pleiotropy.csv/part-00000-d9e99925-4d3d-4d66-b076-a380a977f4ed-c000.csv",data.table=FALSE)
df_betas=open_dataset("data/qsl_l2g_pleiotropy_ta_reasons_of_assoc/")
df_betas=data.frame(df_betas)




unique_study_ids=unique(df_betas$studyId)

n_l2gs=NULL
for (si in unique_study_ids){
  n=sum(df_betas$studyId==si)
  n_l2gs=c(n_l2gs,n)
}
summary(n_l2gs)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    2.00    5.00   15.52   15.00  427.00 

n_threshold=10
studies_to_select=unique_study_ids[n_l2gs>=n_threshold]

out=array(NA,c(length(studies_to_select),6))
out=as.data.frame(out)
colnames(out)=c("studyId","VEP","eQ","pQ","sQ","n")
i=1
for (si in studies_to_select){
  df=df_betas[df_betas$studyId==si,]
  n=nrow(df)
  out[i,1]=si
  out[i,2]=sum(df$VEP==1)/n
  out[i,3]=sum(df$eQTL_coloc==1)/n
  out[i,4]=sum(df$pQTL_coloc==1)/n
  out[i,5]=sum(df$sQTL_coloc==1)/n
  out[i,6]=n
  i=i+1
}
reasons=out

ta_list=c("Haematology","Metabolic","Congenital","Signs.symptoms","Neurology",
          "Immune","Psychiatry","Dermatology",
          "Ophthalmology","Cardiovascular","Oncology","Respiratory","Digestive",
          "Endocrine","Musculoskeletal","Infection","Other")

out=array(NA,c(length(ta_list)+1,10))
out=as.data.frame(out)
colnames(out)=c("ta","VEP","eQ","pQ","sQ","n","VEPp","eQp","pQp","sQp")
i=1
for (ta in ta_list){
  list_ta_studies=unique(df_betas$studyId[df_betas[,ta]==1])
  list_ta_studies=list_ta_studies[list_ta_studies%in%studies_to_select]
  if (length(list_ta_studies)>0){
    df=reasons[(reasons$studyId%in%list_ta_studies),]
    out[i,1]=ta
    out[i,2]=mean(df$VEP)
    out[i,3]=mean(df$eQ)
    out[i,4]=mean(df$pQ)
    out[i,5]=mean(df$sQ)
    out[i,6]=length(list_ta_studies)
    out[i,7]=t.test(df$VEP,reasons$VEP[!(reasons$studyId%in%list_ta_studies)])$p.value
    out[i,8]=t.test(df$eQ,reasons$eQ[!(reasons$studyId%in%list_ta_studies)])$p.value
    out[i,9]=t.test(df$pQ,reasons$pQ[!(reasons$studyId%in%list_ta_studies)])$p.value
    out[i,10]=t.test(df$sQ,reasons$sQ[!(reasons$studyId%in%list_ta_studies)])$p.value
}
  i=i+1
}
out[i,1]="all"
out[i,2]=mean(reasons$VEP)
out[i,3]=mean(reasons$eQ)
out[i,4]=mean(reasons$pQ)
out[i,5]=mean(reasons$sQ)
out[i,6]=nrow(reasons)




