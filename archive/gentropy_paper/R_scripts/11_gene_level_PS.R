setwd("~/Projects/gentropy/notebooks/gentropy_paper/data/")
library(data.table)

gps=fread("20251107_genes_pleiotropy.csv",data.table=FALSE)
gps$year=as.numeric(substr(gps$earliestPublicationDate,1,4))
gps[is.na(gps)]=0
mdl = glm(uniqueDiseases ~ #maxEQTLColoc+maxPQTLColoc+maxVEP+maxDistanceTSS+
            maxEffectiveSampleSize+lofConstraint+misConstraint+geneLength+tissueSpecificity,
          family = quasipoisson(link = "log"),
          data = gps)
summary(mdl)
pred_1 = predict(mdl, type = "response")
cor(pred_1,gps$uniqueDiseases)^2



mdl = glm(uniqueDiseases ~ #maxEQTLColoc+maxPQTLColoc+maxVEP+maxDistanceTSS+
            maxEffectiveSampleSize,
          family = quasipoisson(),
          data = gps)
summary(mdl)
pred_1 = predict(mdl)
cor(pred_1,gps$uniqueDiseases)^2


res=gps$uniqueDiseases-pred_1
res=res-min(res)+1
summary(gps$uniqueDiseases)
summary(res)

mdl = lm(uniqueDiseases ~ maxEffectiveSampleSize,
          data = gps)
pred_1 = predict(mdl)
cor(pred_1,gps$uniqueDiseases)^2





