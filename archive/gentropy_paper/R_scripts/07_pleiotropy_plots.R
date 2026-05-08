library(data.table)

df_v=fread("Projects/gentropy/notebooks/gentropy_paper/data/variant_pleiotropy.csv",data.table=F)
colnames(df_v)

df_g=fread("Projects/gentropy/notebooks/gentropy_paper/data/genes_therapeutic_areas.csv",data.table=F)
colnames(df_g)

l2g_full=fread("Projects/gentropy/notebooks/gentropy_paper/data/l2g_diseases_full.csv",data.table = F)

cs_to_exlude=unique(l2g_full$studyLocusId[duplicated(l2g_full$studyLocusId)])
variants_to_exclude=unique(l2g_full$variantId[l2g_full$studyLocusId%in%cs_to_exlude])
genes_to_exclude=unique(l2g_full$geneId[l2g_full$studyLocusId%in%cs_to_exlude])

df_g=df_g[!(df_g$geneId%in%genes_to_exclude),]
df_v=df_v[!(df_v$variantId%in%variants_to_exclude),]

out=df_g[,c("geneId","uniqueDiseases")]
i=1
for(i in 1:nrow(out)){
  gene=out[i,"geneId"]
  ind=which(l2g_full$geneId==gene)
  snps=unique(l2g_full$variantId[ind])
  ind=which(df_v$variantId%in%snps)
  v_p=df_v$uniqueDiseases[ind]
  out[i,"sum"]=sum(v_p)
  out[i,"max"]=max(v_p)
}


genes_PAV=unique(l2g_full$geneId[l2g_full$VEP==1])
out_PAV=out[out$geneId%in%genes_PAV,]


genes_eQTL=unique(l2g_full$geneId[l2g_full$eQTL_coloc==1])
out_eQTL=out[out$geneId%in%genes_eQTL,]

dat <- data.frame(obs = out[,2], exp1 = out[,3], exp2 = out[,4])
dat <- data.frame(obs = out_PAV[,2], exp1 = out_PAV[,3], exp2 = out_PAV[,4])
dat <- data.frame(obs = out_eQTL[,2], exp1 = out_eQTL[,3], exp2 = out_eQTL[,4])

# Fit models
fit_green <- lm(exp1 ~ obs, data=dat)
fit_red   <- lm(exp2 ~ obs, data=dat)
# X sequence
x_seq <- seq(min(dat$obs), max(dat$obs), length=200)

# Prediction intervals (wider than confidence)
pred_green <- predict(fit_green, newdata=data.frame(obs=x_seq), interval="prediction")
pred_red   <- predict(fit_red,   newdata=data.frame(obs=x_seq), interval="prediction")

par(bty="l", mar=c(5,5,2,2))  # cleaner box

# Plot points
plot(dat$obs, dat$exp1, xlab="Observed pleiotropy", ylab="Expected",
     col=adjustcolor("darkblue", 0.6), pch=16, cex=1.2,xlim=c(0,50),ylim=c(0,150))
points(dat$obs, dat$exp2, col=adjustcolor("darkgreen", 0.6), pch=17, cex=1.2)

# Reference line
abline(a=0, b=1, col="black", lwd=2, lty=3)

# --- GREEN band ---
polygon(c(x_seq, rev(x_seq)),
        c(pred_green[,"lwr"], rev(pred_green[,"upr"])),
        col=adjustcolor("blue", alpha.f=0.15), border=NA)
lines(x_seq, pred_green[,"fit"], col="darkblue", lwd=3)

# --- RED band ---
polygon(c(x_seq, rev(x_seq)),
        c(pred_red[,"lwr"], rev(pred_red[,"upr"])),
        col=adjustcolor("green", alpha.f=0.15), border=NA)
lines(x_seq, pred_red[,"fit"], col="darkgreen", lwd=3)

# Legend
legend("topleft", legend=c("Full independent pleiotropy", "Full shared pleiotropy", "Observed pleiotropy"),
       col=c("darkblue","darkgreen","black"), lwd=c(3,3,2), lty=c(1,1,3),
       pch=c(16,17,NA), pt.cex=1.2, bty="n")













plot(out[,2],out[,3],xlab = "Observed pleiotropy",ylab="Expected",col="green",xlim = c(1,100))
abline(a=lm(out[,3]~out[,2])$coef[1],b=lm(out[,3]~out[,2])$coef[2],col="green")
abline(a=0,b=1,col="black")
points(out[,2],out[,4],col="red")
abline(a=lm(out[,4]~out[,2])$coef[1],b=lm(out[,4]~out[,2])$coef[2],col="red")






# Rename for clarity
dat <- data.frame(obs = out[,2], exp1 = out[,3], exp2 = out[,4])

# Fit models
fit_green <- lm(exp1 ~ obs, data=dat)
fit_red   <- lm(exp2 ~ obs, data=dat)

# X sequence
x_seq <- seq(min(dat$obs), max(dat$obs), length=200)

# Prediction intervals (wider than confidence)
pred_green <- predict(fit_green, newdata=data.frame(obs=x_seq), interval="prediction")
pred_red   <- predict(fit_red,   newdata=data.frame(obs=x_seq), interval="prediction")

# Plot
plot(dat$obs, dat$exp1, xlab="Observed pleiotropy", ylab="Expected",
     col="green", xlim=c(1,100))
points(dat$obs, dat$exp2, col="red")

# Reference line
abline(a=0, b=1, col="black")

# --- GREEN band ---
polygon(c(x_seq, rev(x_seq)),
        c(pred_green[,"lwr"], rev(pred_green[,"upr"])),
        col=adjustcolor("green", alpha.f=0.2), border=NA)
lines(x_seq, pred_green[,"fit"], col="green", lwd=2)

# --- RED band ---
polygon(c(x_seq, rev(x_seq)),
        c(pred_red[,"lwr"], rev(pred_red[,"upr"])),
        col=adjustcolor("red", alpha.f=0.2), border=NA)
lines(x_seq, pred_red[,"fit"], col="red", lwd=2)


######

# Basic style


hist(df$uniqueDiseases)

plot(df$maxMAF,df$uniqueDiseases)

bins=c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
out=array(NA,c(length(bins),3))


i=1
for (i in 1:length(bins)){
  thr2=bins[i]
  if (i!=1) thr1=bins[i-1] else thr1=0
  ind=which(df$maxMAF>=thr1 & df$maxMAF<thr2)
  out[i,1]=thr2
  out[i,2]=mean(df$uniqueDiseases[ind])
  out[i,3]=sd(df$uniqueDiseases[ind])/sqrt(length(ind))
}


summary(lm(df$uniqueDiseases~df$maxMAF))

plot(df$maxMAF,df$uniqueDiseases,ylim=c(1,2))
l=lm(df$uniqueDiseases~df$maxMAF+I(df$maxMAF^2))
summary(l)
points(df$maxMAF,predict(l),col='green')


df1=df[df$vepScore<0.66,]
l=lm(df1$uniqueDiseases~df1$maxMAF+I(df1$maxMAF^2))
summary(l)
points(df1$maxMAF,predict(l),col='red')


l=lm(df$uniqueDiseases~df$maxAbsBeta+I(df$maxAbsBeta^2))
summary(l)
plot(df$maxAbsBeta,predict(l))


l=lm(df$uniqueDiseases~df$maxAbsBeta*df$maxMAF*I(df$vepScore>=0.66))
summary(l)




