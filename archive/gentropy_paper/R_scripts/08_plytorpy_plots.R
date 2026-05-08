setwd("~/Desktop/PSs/")
library(data.table)

vlp=fread("variants_with_fitness_scores_pred.csv",data.table=F)
vlp$weighted_sum_beta_variant[is.na(vlp$weighted_sum_beta_variant)]=0.0001
vlp$weighted_sum_beta_variant=abs(vlp$weighted_sum_beta_variant)

vlp_pl=vlp[vlp$uniqueDiseases>1,]

vlp$res1=pchisq(32.84125,df=1,low=F,ncp=(vlp$maxAbsBeta^2)*vlp$maxEffectiveSampleSize*vlp$maxVarG/11)
vlp$res1=(vlp$res1/mean(vlp$res1))*mean(vlp$uniqueDiseases)





library(MASS)

mdl = glm(uniqueDiseases ~ I(maxMAF *(maxAbsBeta^2) * maxEffectiveSampleSize),
          family = quasipoisson(link = "log"),
          data = vlp)
vlp$stat_pred = predict(mdl, type = "response")
cor(vlp$stat_pred,vlp$uniqueDiseases)^2

mdl = glm(uniqueDiseases ~ res1,
          family = quasipoisson(link = "log"),
          data = vlp)
vlp$stat_res = predict(mdl, type = "response")
cor(vlp$stat_res,vlp$uniqueDiseases)^2

mdl = glm(uniqueDiseases ~ res1+I(maxMAF *(maxAbsBeta^2) * maxEffectiveSampleSize),
          family = quasipoisson(link = "log"),
          data = vlp)
vlp$stat_res = predict(mdl, type = "response")
cor(vlp$stat_res,vlp$uniqueDiseases)^2


mdl = glm(uniqueDiseases ~ res1,
          family = quasipoisson(link = "log"),
          data = vlp)
vlp$stat_res = mdl$residuals
vlp$stat_res=(vlp$stat_res/mean(vlp$stat_res))*mean(vlp$uniqueDiseases)
cor(vlp$stat_res,vlp$uniqueDiseases)

mdl = glm(uniqueDiseases ~ maxEffectiveSampleSize,
          family = quasipoisson(link = "log"),
          data = vlp)
vlp$stat_res = mdl$residuals
vlp$stat_res=(vlp$stat_res/mean(vlp$stat_res))*mean(vlp$uniqueDiseases)
cor(vlp$stat_res,vlp$uniqueDiseases)

#######
library(dplyr)
library(ggplot2)

ind_pleiotropy=which(vlp$uniqueDiseases>1) 


# Define the bin breaks
maf_breaks <- c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)  # 0, 0.01, 0.02, …, 0.5
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp=vlp[vlp$uniqueDiseases>=1,]

vlp_summary <-  tmp %>%
  mutate(
    maxMAF_bin = cut(maxMAF,
                     breaks = maf_breaks,
                     include.lowest = TRUE,
                     labels = maf_labels)
  ) %>%
  group_by(maxMAF_bin) %>%
  summarise(
    mean_uniqueDiseases = mean(res1, na.rm = TRUE),
    se_uniqueDiseases   = sd(res1, na.rm = TRUE) / sqrt(sum(!is.na(res1))),
    lower = mean_uniqueDiseases - 1.96 * se_uniqueDiseases,
    upper = mean_uniqueDiseases + 1.96 * se_uniqueDiseases,
    .groups = "drop"
  )

# Plot
ggplot(vlp_summary, aes(x = maxMAF_bin, y = mean_uniqueDiseases)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.002) +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Number of traits",
    title = "Modeling of pleytopy "
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




####
tmp=vlp
library(dplyr)
library(ggplot2)
library(tidyr)

# Define MAF bins and labels
maf_breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

# Filter data
tmp <- vlp %>% dplyr::filter(uniqueDiseases >= 1)

# Summarize only the two variables of interest
vlp_summary <- tmp %>%
  dplyr::mutate(
    maxMAF_bin = cut(maxMAF,
                     breaks = maf_breaks,
                     include.lowest = TRUE,
                     labels = maf_labels)
  ) %>%
  dplyr::group_by(maxMAF_bin) %>%
  dplyr::summarise(
    mean_uniqueDiseases = mean(uniqueDiseases, na.rm = TRUE),
    mean_res1           = mean(res1, na.rm = TRUE),
    se_uniqueDiseases   = sd(uniqueDiseases, na.rm = TRUE) / sqrt(sum(!is.na(uniqueDiseases))),
    se_res1             = sd(res1, na.rm = TRUE) / sqrt(sum(!is.na(res1))),
    .groups = "drop"
  )

# Convert to long format (only two lines)
vlp_long <- vlp_summary %>%
  tidyr::pivot_longer(
    cols = c(mean_uniqueDiseases, mean_res1),
    names_to = "variable",
    values_to = "mean_value"
  ) %>%
  dplyr::mutate(
    se = dplyr::case_when(
      variable == "mean_uniqueDiseases" ~ se_uniqueDiseases,
      variable == "mean_res1" ~ se_res1
    ),
    lower = mean_value - 1.96 * se,
    upper = mean_value + 1.96 * se,
    variable = dplyr::recode(variable,
                             mean_uniqueDiseases = "Observed (uniqueDiseases)",
                             mean_res1 = "Predicted")
  )

# Plot the two lines
ggplot(vlp_long, aes(x = maxMAF_bin, y = mean_value, color = variable, group = variable)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02, alpha = 0.6) +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Number of traits",
    title = "Observed vs Predicted Number of Traits by MAF Bin",
    color = "Metric"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#####

























library(dplyr)
library(ggplot2)

# Define bins for uniqueDiseases (1–10+)
disease_breaks <- c(0, 1, 2, 3, 4,5, Inf)
disease_labels <- c("1", "2", "3", "4","5", "5+")

# Filter if needed
tmp <- vlp[vlp$uniqueDiseases >= 0, ]

# Summarise mean maxAbsBeta by uniqueDiseases bins
vlp_summary <- tmp %>%
  mutate(
    uniqueDiseases_bin = cut(uniqueDiseases,
                             breaks = disease_breaks,
                             include.lowest = TRUE,
                             right = FALSE,
                             labels = disease_labels)
  ) %>%
  group_by(uniqueDiseases_bin) %>%
  summarise(
    mean_maxAbsBeta = mean(maxAbsBeta, na.rm = TRUE),
    se_maxAbsBeta   = sd(maxAbsBeta, na.rm = TRUE) / sqrt(sum(!is.na(maxAbsBeta))),
    lower = mean_maxAbsBeta - 1.96 * se_maxAbsBeta,
    upper = mean_maxAbsBeta + 1.96 * se_maxAbsBeta,
    .groups = "drop"
  )

# Plot mean maxAbsBeta ± 95% CI by uniqueDiseases_bin
ggplot(vlp_summary, aes(x = uniqueDiseases_bin, y = mean_maxAbsBeta)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_minimal() +
  labs(
    x = "Unique Diseases (binned)",
    y = "Mean |β| ± 95% CI",
    title = "Mean maxAbsBeta by Unique Diseases bin"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



library(dplyr)
library(tidyr)
library(ggplot2)

maf_breaks <- c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp <- vlp[vlp$uniqueDiseases >= 1, ]

vlp_summary <- tmp %>%
  mutate(
    maxMAF_bin = cut(maxMAF, breaks = maf_breaks, include.lowest = TRUE, labels = maf_labels)
  ) %>%
  dplyr::select(maxMAF_bin, uniqueDiseases, stat_pred) %>%
  tidyr::pivot_longer(
    cols = c(uniqueDiseases, stat_pred),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::group_by(maxMAF_bin, variable) %>%
  dplyr::summarise(
    mean_val = mean(value, na.rm = TRUE),
    se_val   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    lower = mean_val - 1.96 * se_val,
    upper = mean_val + 1.96 * se_val,
    .groups = "drop"
  )



ggplot(vlp_summary, aes(x = maxMAF_bin, y = mean_val, color = variable, group = variable)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.002) +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Mean ± 95% CI",
    title = "Mean observed vs predicted by maxMAF bin",
    color = "Variable"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



