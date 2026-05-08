setwd("~/Projects/gentropy/notebooks/gentropy_paper/")

library(data.table)
library(arrow)
x=fread("data/protein_coding_eqtls_l2g_fm_based_rare_varint_pleiotropy.csv/part-00000-d9e99925-4d3d-4d66-b076-a380a977f4ed-c000.csv",data.table=FALSE)
#df_betas=open_dataset("data/qsl_l2g_pleiotropy/")
#df_betas=data.frame(df_betas)

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

x$VEP_low_pl=0
x$VEP_low_pl[x$VEP==1 & x$pleiotropy<0.2]=1

x$VEP_high_pl=0
x$VEP_high_pl[x$VEP==1 & x$pleiotropy>=0.2]=1

x$pqtl_high_pl=0
x$pqtl_high_pl[x$pleiotropy>=0.2 & x$pQTL_coloc==1]=1

x$pqtl_low_pl=0
x$pqtl_low_pl[x$pleiotropy<0.2 & x$pQTL_coloc==1]=1
####
####
####
colnames_to_impute=c("eQTL_coloc","l2g","VEP",
                     "pQTL_coloc","sQTL_coloc","distance","rare_variant_evidence","VEP_low_pl","VEP_high_pl","pqtl_high_pl","pqtl_low_pl")
colnames_to_comapre=c(colnames_to_impute,"maxClinicalTrialPhase","hasSafetyEvent")

out=array(NA,c(length(colnames_to_impute)+4,3))
out=as.data.frame(out)
c1="geneticConstraint"
c2="tissueDistribution"
colnames(out)=c("Category",c1,c2)
l=1
for (i in colnames_to_impute){
  out[l,1]=i
  out[l,2]=mean(x[,c1][x[,i]==1],na.rm=T)
  out[l,3]=mean(x[,c2][x[,i]==1],na.rm=T)
  l=l+1
}
out[l,1]="all"
out[l,2]=mean(x[,c1],na.rm=T)
out[l,3]=mean(x[,c2],na.rm=T)
l=l+1


out[l,1]="chembl"
out[l,2]=mean(x[,c1][x$maxClinicalTrialPhase>=0.5],na.rm=T)
out[l,3]=mean(x[,c2][x$maxClinicalTrialPhase>=0.5],na.rm=T)
l=l+1

out[l,1]="safety"
out[l,2]=mean(x[,c1][x$hasSafetyEvent==1],na.rm=T)
out[l,3]=mean(x[,c2][x$hasSafetyEvent==1],na.rm=T)
l=l+1

out[l,1]="eQTLs"
out[l,2]=mean(x[,c1][x$eqtl==1],na.rm=T)
out[l,3]=mean(x[,c2][x$eqtl==1],na.rm=T)
l=l+1


library(ggplot2)

data=out
colnames(data)=c("Category","c1","c2")

# Basic scatter plot with text labels
ggplot(data, aes(x = c1, y = c2, label = Category)) +
  geom_point() +          # plot points
  geom_text(vjust = -0.5, hjust = 0.5) +  # add labels slightly above points
  theme_minimal() +       # cleaner theme
  labs(x = c1, y = c2, 
       title = "")



####
####
####
colnames_to_impute=c("eQTL_coloc","l2g","VEP",
                     "pQTL_coloc","sQTL_coloc","distance","rare_variant_evidence","VEP_low_pl","VEP_high_pl","pqtl_high_pl","pqtl_low_pl")
colnames_to_comapre=c(colnames_to_impute,"maxClinicalTrialPhase","hasSafetyEvent")

out=array(NA,c(length(colnames_to_impute),3))
out=as.data.frame(out)
c1="overlap with chembl"
c2="tissueDistribution"
colnames(out)=c("Category",c1,c2)
l=1

chembl_genes=unique(x$geneId[x$maxClinicalTrialPhase>=0.5])

for (i in colnames_to_impute){
  out[l,1]=i
  out[l,2]=sum(x$geneId[x[,i]==1] %in% chembl_genes)/length(x$geneId[x[,i]==1])
  out[l,3]=mean(x[,c2][x[,i]==1],na.rm=T)
  l=l+1
}


library(ggplot2)

data=out
colnames(data)=c("Category","c1","c2")

# Basic scatter plot with text labels
ggplot(data, aes(x = c1, y = c2, label = Category)) +
  geom_point() +          # plot points
  geom_text(vjust = -0.5, hjust = 0.5) +  # add labels slightly above points
  theme_minimal() +       # cleaner theme
  labs(x = c1, y = c2, 
       title = "")





####
####
####
colnames_to_impute=c("eQTL_coloc","l2g","VEP",
                     "pQTL_coloc","sQTL_coloc","distance","rare_variant_evidence")
colnames_to_comapre=c(colnames_to_impute,"maxClinicalTrialPhase","hasSafetyEvent")

out=array(NA,c(length(colnames_to_impute)+4,3))
out=as.data.frame(out)
c1="geneticConstraint"
c2="distinct_biosample_count"
colnames(out)=c("Category",c1,c2)
l=1
for (i in colnames_to_impute){
  out[l,1]=i
  out[l,2]=mean(x[,c1][x[,i]==1],na.rm=T)
  out[l,3]=mean(x[,c2][x[,i]==1],na.rm=T)
  l=l+1
}
out[l,1]="all"
out[l,2]=mean(x[,c1],na.rm=T)
out[l,3]=mean(x[,c2],na.rm=T)
l=l+1


out[l,1]="chembl"
out[l,2]=mean(x[,c1][x$maxClinicalTrialPhase>=0.5],na.rm=T)
out[l,3]=mean(x[,c2][x$maxClinicalTrialPhase>=0.5],na.rm=T)
l=l+1

out[l,1]="safety"
out[l,2]=mean(x[,c1][x$hasSafetyEvent==1],na.rm=T)
out[l,3]=mean(x[,c2][x$hasSafetyEvent==1],na.rm=T)
l=l+1

out[l,1]="eQTLs"
out[l,2]=mean(x[,c1][x$eqtl==1],na.rm=T)
out[l,3]=mean(x[,c2][x$eqtl==1],na.rm=T)
l=l+1


library(ggplot2)

data=out
colnames(data)=c("Category","c1","c2")

# Basic scatter plot with text labels
ggplot(data, aes(x = c1, y = c2, label = Category)) +
  geom_point() +          # plot points
  geom_text(vjust = -0.5, hjust = 0.5) +  # add labels slightly above points
  theme_minimal() +       # cleaner theme
  labs(x = c1, y = c2, 
       title = "")





####
#### PLEYOTROPY VS TISSUE DISTRIBUTION
####
colnames_to_impute=c("eQTL_coloc","l2g","VEP",
                     "pQTL_coloc","sQTL_coloc","distance")
out=array(NA,c(length(colnames_to_impute),3))
out=as.data.frame(out)
c1="pleiotropy"
c2="tissueDistribution"
colnames(out)=c("Category",c1,c2)
l=1
xx=x[x$pleiotropy>0,]
for (i in colnames_to_impute){
  out[l,1]=i
  out[l,2]=mean(xx[,c1][xx[,i]==1],na.rm=T)
  out[l,3]=mean(xx[,c2][xx[,i]==1],na.rm=T)
  l=l+1
}


library(ggplot2)

data=out
colnames(data)=c("Category","c1","c2")

# Basic scatter plot with text labels
ggplot(data, aes(x = c1, y = c2, label = Category)) +
  geom_point() +          # plot points
  geom_text(vjust = -0.5, hjust = 0.5) +  # add labels slightly above points
  theme_minimal() +       # cleaner theme
  labs(x = c1, y = c2, 
       title = "")










#####
##### IGNORE BELOW
#####


colnames(x)


summary(lm(distinct_biosample_count~tissueSpecificity+tissueDistribution+geneticConstraint,x))

summary(lm(geneticConstraint~l2g,x))

table(x$l2g==1 & x$distinct_biosample_count>0)

y=x[x$l2g==1,]
t.test(y$geneticConstraint[y$eqtl==0],y$geneticConstraint[y$eqtl==1])

y=x[x$eqtl==1,]



summary(lm(y$geneticConstraint~y$tissueDistribution+y$tissueSpecificity))
summary(lm(y$geneticConstraint~y$distinct_biosample_count))
summary(lm(x$hasSafetyEvent~x$distinct_biosample_count+x$tissueDistribution+x$tissueSpecificity+x$geneticConstraint+x$l2g))




summary(lm(geneticConstraint~distinct_biosample_count+VEP+pQTL_coloc+eQTL_coloc+distance+tissueSpecificity+rare_variant_evidence+maxClinicalTrialPhase+hasSafetyEvent+sQTL_coloc+pleiotropy,data=x))



summary(lm(hasSafetyEvent~distinct_biosample_count+geneticConstraint+VEP+pQTL_coloc+eQTL_coloc+distance+tissueSpecificity+rare_variant_evidence+l2g+sQTL_coloc,data=x))


summary(lm(hasSafetyEvent~geneticConstraint+pQTL_coloc+eQTL_coloc+distance+tissueSpecificity,data=x))


summary(lm(geneticConstraint~distinct_biosample_count+pQTL_coloc+eQTL_coloc+distance+tissueSpecificity+rare_variant_evidence+maxClinicalTrialPhase+hasSafetyEvent+pleiotropy+VEP,data=x))

summary(lm(pleiotropy~distinct_biosample_count+pQTL_coloc+eQTL_coloc+distance+tissueSpecificity+VEP,data=x[x$pleiotropy>0,]))

summary(lm(geneticConstraint~tissueSpecificity+pleiotropy*VEP,data=x[x$l2g>0,]))

summary(lm(hasSafetyEvent~pleiotropy,data=x))

summary(lm(pleiotropy~VEP+eQTL_coloc+pQTL_coloc+sQTL_coloc+distinct_biosample_count,data=x[x$l2g>0,]))







