library(data.table)

pleiotropy=fread("Projects/gentropy/notebooks/gentropy_paper/data/pleiotropy_combined_evidence.csv",data.table=F)
combined_evidence=fread("Projects/gentropy/notebooks/gentropy_paper/data/combined_evidence.csv",data.table = F)
target=fread("Projects/gentropy/notebooks/gentropy_paper/data/target_with_constraints_2509.csv",data.table = F)

target$lof_constr[is.na(target$lof_constr)]=mean(target$lof_constr,na.rm=T)
#target$lof_constr=scale(target$lof_constr)
target$lof_constr=-(target$lof_constr)

table(combined_evidence$datasourceId)

length(unique(combined_evidence$targetId[combined_evidence$datasourceId=="l2g_combined"]))
length(unique(combined_evidence$targetId[combined_evidence$datasourceId=="gene_burden"]))
length(unique(combined_evidence$targetId[combined_evidence$datasourceId=="rare"]))

rare=unique(combined_evidence$targetId[combined_evidence$datasourceId=="rare"])
gene_burden=unique(combined_evidence$targetId[combined_evidence$datasourceId=="gene_burden"])
l2g_combined=unique(combined_evidence$targetId[combined_evidence$datasourceId=="l2g_combined"])
l2g_PAV=unique(combined_evidence$targetId[combined_evidence$datasourceId=="l2g_VEP"])
l2g_eQTL=unique(combined_evidence$targetId[combined_evidence$datasourceId=="l2g_eQTL"])
inter=unique(c(intersect(rare,gene_burden),intersect(rare,l2g_combined),intersect(gene_burden,l2g_combined)))
only_rare=rare[!(rare%in%inter)]
only_gene_burden=gene_burden[!(gene_burden%in%inter)]
only_l2g_combined=l2g_combined[!(l2g_combined%in%inter)]

l2g_PAV=l2g_PAV[!(l2g_PAV%in%inter)]
l2g_eQTL=l2g_eQTL[!(l2g_eQTL%in%inter)]
l2g_eQTL=l2g_eQTL[!(l2g_eQTL%in%l2g_PAV)]



to_plot=as.data.frame(
      array(NA,
        c(
          length(rare)+
            length(gene_burden)+
            length(l2g_combined)+
            length(inter)+
            length(only_gene_burden)+
            length(only_rare)+
            length(only_l2g_combined)+
            length(l2g_PAV)+
            length(l2g_eQTL),
          4
          )
        )
)
colnames(to_plot)=c("dataSource","pleiotropy","constrain","geneId")
constr="lof_constr"

df=rare
ttl="Rare diseases"
ind_to_plot_1=1
ind_to_plot_2=length(df)
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=l2g_combined
ttl="GWAS"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=gene_burden
ttl="Gene Burden"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=inter
ttl="Intrsect"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=only_gene_burden
ttl="only_gene_burden"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=only_rare
ttl="only_rare"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=only_l2g_combined
ttl="only_l2g_combined"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=l2g_PAV
ttl="l2g_PAV"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

df=l2g_eQTL
ttl="l2g_eQTL"
ind_to_plot_1=ind_to_plot_2+1
ind_to_plot_2=ind_to_plot_1+length(df)-1
ind_pl=match(df,pleiotropy$targetId)
ind_tr=match(df,target$targetId)

to_plot[ind_to_plot_1:ind_to_plot_2,1]=ttl
to_plot[ind_to_plot_1:ind_to_plot_2,4]=df
to_plot[ind_to_plot_1:ind_to_plot_2,2]=pleiotropy$unique_disease_count[ind_pl]
to_plot[ind_to_plot_1:ind_to_plot_2,3]=target[ind_tr,constr]

library(dplyr)
library(ggplot2)

# Summarise: mean, SE, and 95% CI for each group
summary_df <- to_plot %>%
  group_by(dataSource) %>%
  summarise(
    mean_pleiotropy = mean(pleiotropy, na.rm = TRUE),
    se_pleiotropy = sd(pleiotropy, na.rm = TRUE)/sqrt(n()),
    mean_constrain = mean(constrain, na.rm = TRUE),
    se_constrain = sd(constrain, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ci_pleiotropy = 1.96 * se_pleiotropy,
    ci_constrain = 1.96 * se_constrain
  )

print(summary_df)

summary_df_small=summary_df[-c(1,2,4),]

# Plot: x = mean_pleiotropy, y = mean_constrain
ggplot(summary_df_small, aes(x = mean_pleiotropy, y = mean_constrain, color = dataSource)) +
  geom_point(size = 3) +
  # error bars for pleiotropy (x-axis)
  geom_errorbarh(aes(xmin = mean_pleiotropy - ci_pleiotropy, 
                     xmax = mean_pleiotropy + ci_pleiotropy), height = 0.005) +
  # error bars for constrain (y-axis)
  geom_errorbar(aes(ymin = mean_constrain - ci_constrain, 
                    ymax = mean_constrain + ci_constrain), width = 0.05) +
  theme_minimal() +
  labs(x = "Mean pleiotropy (±95% CI)", 
       y = "Mean constraint (±95% CI)")




library(ggplot2)

# Assuming your data is in a data frame called to_plot
ggplot(to_plot, aes(x = dataSource, y = pleiotropy, fill = dataSource)) +
  geom_violin(trim = FALSE, alpha = 0.6) +   # violin shape
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +  # boxplot inside
  geom_jitter(width = 0.1, alpha = 0.01, size = 2) +  # raw points
  theme_minimal() +
  ylim(0, 10)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    x = "Data Source",
    y = "Pleiotropy",
    title = "Distribution of Pleiotropy by Data Source"
  )


#####
length(to_plot$geneId[to_plot$dataSource=="Intrsect"])
cg=fread("Projects/gentropy/notebooks/gentropy_paper/data/genes_with_10_or_more_diseases.txt",data.table=F,header = F)[,1]
length(cg)
table(cg%in%to_plot$geneId[to_plot$dataSource=="Intrsect"])
