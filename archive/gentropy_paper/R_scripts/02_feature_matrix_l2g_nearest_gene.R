setwd("~/Projects/gentropy/notebooks/gentropy_paper/")

library(data.table)
x=fread("data/feature_matrix_l2g.csv/part-00000-b83911c7-398c-4144-8ca6-25fa0d48674c-c000.csv",data.table=FALSE)
colnames(x)
#[1] "geneId"                                 "studyLocusId"                           "credibleSetConfidence"                 
#[4] "distanceFootprintMean"                  "distanceFootprintMeanNeighbourhood"     "distanceSentinelFootprint"             
#[7] "distanceSentinelFootprintNeighbourhood" "distanceSentinelTss"                    "distanceSentinelTssNeighbourhood"      
#[10] "distanceTssMean"                        "distanceTssMeanNeighbourhood"           "eQtlColocClppMaximum"                  
#[13] "eQtlColocClppMaximumNeighbourhood"      "eQtlColocH4Maximum"                     "eQtlColocH4MaximumNeighbourhood"       
#[16] "geneCount500kb"                         "isProteinCoding"                        "pQtlColocClppMaximum"                  
#[19] "pQtlColocClppMaximumNeighbourhood"      "pQtlColocH4Maximum"                     "pQtlColocH4MaximumNeighbourhood"       
#[22] "proteinGeneCount500kb"                  "sQtlColocClppMaximum"                   "sQtlColocClppMaximumNeighbourhood"     
#[25] "sQtlColocH4Maximum"                     "sQtlColocH4MaximumNeighbourhood"        "vepMaximum"                            
#[28] "vepMaximumNeighbourhood"                "vepMean"                                "vepMeanNeighbourhood"                  
#[31] "score"                                  "l2g_05"                                 "eQTL_coloc"                            
#[34] "pQTL_coloc"                             "sQTL_coloc"                             "VEP"                                   
#[37] "distance"   

y=x[x$score>=0.5,]
table(y$distance==0 & y$eQTL_coloc==1)


table(x$distance==1 & x$eQTL_coloc==1)

l2g05=unique(y$studyLocusId)
all_slids=unique(x$studyLocusId)
no_05=all_slids[!all_slids%in%l2g05]

lll=rep(0,nrow(x))
for(i in all_slids){
  ind=which(x$studyLocusId==i)
  lll[ind[which.max(x$score[ind])[1]]]=1
}


# Convert to data.table
dt <- as.data.table(x)

# Initialize lll
lll <- rep(0, nrow(dt))

# For each studyLocusId, find the row with max score
max_rows <- dt[, .I[which.max(score)], by = studyLocusId]$V1

# Set lll
lll[max_rows] <- 1

y=x[lll==1 | x$score>=0.5,]
table(y$distance)


tra_set=fread("../2503_release/training_set/patched_training_2503-testrun-1_all_string_005_extended_EGL_variants.tsv",data.table=FALSE)
posit=tra_set[tra_set$GSP==1,]

keys_x=paste(x$geneId,x$studyLocusId,sep = "_")
keys_pos=paste(posit$geneId,posit$studyLocusId,sep = "_")

table(keys_pos%in%keys_x)
ind=match(keys_pos,keys_x)
table(keys_x[ind]==keys_pos)

y=x[ind,]



