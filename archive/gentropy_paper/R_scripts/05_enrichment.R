library(data.table)

df=fread("Projects/gentropy/notebooks/gentropy_paper/data/data_for_drug_enrichment_ta.csv",data.table=F)
colnames(df)
#[1] "targetId"                   "diseaseId"                  "indirect_assoc_score"       "max_nSamples"              
#[5] "yearsSinceFirstPublication" "eQTL_coloc"                 "pQTL_coloc"                 "VEP"                       
#[9] "distanceTSS"                "maxClinicalPhase"           "geneticSupport"             "outcome"   

df=df[df$totalTherapeuticAreas!=0,]
dim(df)

table(df$mappedTherapeuticAreas)

df$maxNSamples_scaled <- scale(df$maxNSamples)
df$mappedTherapeuticAreas=as.factor(df$mappedTherapeuticAreas)

L0=summary(glm(outcome~geneticSupport,data=df,family="binomial"))
L0

L2=summary(glm(outcome~geneticSupport+maxNSamples_scaled,data=df,family="binomial"))
L2

library(lme4)
model <- glmer(outcome~geneticSupport+(1|mappedTherapeuticAreas),
               data=df,
               family = binomial(link = "logit"))
summary(model)

m1 <- glmer(outcome ~ geneticSupport +maxNSamples_scaled+ (1 | mappedTherapeuticAreas),
            data = df, family = binomial)
summary(m1)




m1 <- glmer(outcome ~ geneticSupport + (1 | mappedTherapeuticAreas),
            data = df, family = binomial)

m0 <- glm(outcome ~ geneticSupport, 
            data = df, family = binomial)

lrt <- 2 * (logLik(m1) - logLik(m0))
p <- 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE)





df$max_nSamples_scaled <- scale(df$max_nSamples)
df$yearsSinceFirstPublication_scaled <- scale(df$yearsSinceFirstPublication)

m1 <- glmer(outcome ~ geneticSupport +
              max_nSamples_scaled +
              yearsSinceFirstPublication_scaled +
              (1 | mappedTherapeuticAreas),
            data = df, family = binomial)
summary(m1)


m1 <- glm(outcome ~ geneticSupport +
              max_nSamples_scaled +
              yearsSinceFirstPublication_scaled
            data = df, family = binomial)
summary(m1)


m1 <- glmer(outcome ~ geneticSupport +
              (1 | mappedTherapeuticAreas),
            data = df, family = binomial)
summary(m1)
