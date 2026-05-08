setwd("~/Projects/gentropy/notebooks/gentropy_paper/")

library(data.table)
library(arrow)
x=fread("data/protein_coding_eqtls_l2g_fm_based_rare_varint_pleiotropy.csv/part-00000-d9e99925-4d3d-4d66-b076-a380a977f4ed-c000.csv",data.table=FALSE)
df_betas=open_dataset("data/qsl_l2g_pleiotropy/")
df_betas=data.frame(df_betas)


######
table(x$hasSafetyEvent)
table(x$maxClinicalTrialPhase)

x$distinct_biosample_count[is.na(x$distinct_biosample_count)]=0
x$l2g[is.na(x$l2g)]=0

x$hasSafetyEvent[is.na(x$hasSafetyEvent)]=0
x$maxClinicalTrialPhase[is.na(x$maxClinicalTrialPhase)]=0


colnames_to_impute=c("eQTL_coloc","l2g","VEP",
                     "pQTL_coloc","sQTL_coloc","distance","rare_variant_evidence","geneticConstraint","pleiotropy")
for (i in colnames_to_impute){
  x[is.na(x[,i]),i]=0
}

colnames_to_zero=c("eQTL_coloc","VEP",
                   "pQTL_coloc","sQTL_coloc","distance")
for (i in colnames_to_zero){
  x[x$l2g==0,i]=0
}


x$eqtl=0
x$eqtl[x$distinct_biosample_count>0]=1

x$hasSafetyEvent[x$hasSafetyEvent==-1]=1
x$tissueDistribution=-x$tissueDistribution
x$geneticConstraint=-x$geneticConstraint

####

df_betas$lofConstraint=-df_betas$lofConstraint
df_betas$tissueDistribution=-df_betas$tissueDistribution

colnames_to_include=c("eQTL_coloc","VEP",
                   "pQTL_coloc","sQTL_coloc","distance")

for (cln in colnames_to_include){
  df_betas[,cln]=0
  list_of_genes=unique(x$geneId[x[,cln]==1])
  df_betas[df_betas$geneId%in%list_of_genes,cln]=1
}

df_betas[,"chembl"]=0
list_of_genes=unique(x$geneId[x$maxClinicalTrialPhase>=0.5])
df_betas[df_betas$geneId%in%list_of_genes,"chembl"]=1

df_betas[,"safety"]=0
list_of_genes=unique(x$geneId[x$hasSafetyEvent==1])
df_betas[df_betas$geneId%in%list_of_genes,"safety"]=1

df_betas[,"abs_beta"]=abs(df_betas$rescaledStatistics$estimatedBeta)
df_betas[,"prev"]=abs(df_betas$rescaledStatistics$prev)
df_betas[,"chi2"]=abs(df_betas$rescaledStatistics$chi2Stat)
df_betas[,"vep_score"]=abs(df_betas$vepEffect$normalisedScore)

columns_for_pca=c("nSamples","vep_score",
                  "majorPopulationMAF","abs_beta","prev","chi2",
                  "lofConstraint","misConstraint","synConstraint","pleiotropy","tissueDistribution",
                  "eQTL_coloc","VEP","pQTL_coloc","sQTL_coloc","distance",
                  "chembl","safety","geneId")
df=df_betas[,columns_for_pca]

#dim(df)
#dim(na.omit(df))
#pca_result <- prcomp(na.omit(df), center = TRUE, scale. = TRUE)
#plot(pca_result, type = "l")

#biplot(pca_result)


data=na.omit(df)
data <- unique(data)
gene_id_labels=data$geneId
chembl=data$chembl
safety=data$safety
data=data[,-which(colnames(data)%in%c("geneId","chembl","safety"))]
data=scale(data)
data=as.data.table(data)
#summary(lm(chembl~data$eQTL_coloc+data$VEP*data$pleiotropy+data$abs_beta+data$pQTL_coloc+data$lofConstraint+data$abs_beta))
#summary(lm(chembl~as.matrix(data)+data$VEP*data$pleiotropy))
#summary(lm(safety~data))

summary(lm(chembl~data$majorPopulationMAF+data$VEP*data$pleiotropy*data$pQTL_coloc+data$eQTL_coloc*data$tissueDistribution+data$lofConstraint+data$misConstraint
           +data$synConstraint+data$sQTL_coloc+data$distance))

summary(lm(chembl~data$VEP*data$pleiotropy))

#######

ind=which(chembl==1)
summary(lm(safety[ind]~data$majorPopulationMAF[ind]+data$VEP[ind]*data$pleiotropy[ind]*data$pQTL_coloc[ind]+data$eQTL_coloc[ind]*data$tissueDistribution[ind]))



library(umap)
umap_result <- umap(data,)
chembl_col=chembl
chembl_col[chembl==0]=10
plot(umap_result$layout, pch = 1, main = "UMAP projection",col=chembl_col)

library(Rtsne)
library(ggplot2)


set.seed(42)  # For reproducibility
tsne_result <- Rtsne(data, dims = 2, perplexity = 5000, theta = 0.9, verbose = TRUE, max_iter = 500)

tsne_df <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  Label = gene_id_labels
)
plot(tsne_df[,1:2])



