setwd("~/Desktop/PSs/")
library(data.table)

vlp=fread("variants_with_fitness_scores.csv",data.table=F)

library(MASS)

# Example model
vlp <- vlp %>%
  mutate(
    pred_var = log1p(maxEffectiveSampleSize*maxVarG*maxAbsBeta^2) # log(1+x) avoids log(0)
  )

fit_nb <- glm.nb(uniqueDiseases ~ I(pred_var+1), data = vlp)
hist(predict(fit_nb))
cor(predict(fit_nb),vlp$uniqueDiseases)

vlp$res1=vlp$uniqueDiseases-predict(fit_nb)


library(dplyr)
library(ggplot2)

# Define MAF bins and labels
maf_breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp <- vlp %>%
  mutate(
    maxMAF_bin = cut(
      maxMAF,
      breaks = maf_breaks,
      include.lowest = TRUE,
      labels = maf_labels
    )
  )

# Count number of variants per bin
maf_counts <- tmp %>%
  group_by(maxMAF_bin) %>%
  summarise(
    N_variants = n(),
    .groups = "drop"
  )

# Bar plot
ggplot(maf_counts, aes(x = maxMAF_bin, y = N_variants)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Number of variants",
    title = "Number of variants per MAF bin"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#Pleytopicity - proportions per maf bins

# Define MAF bins and labels
maf_breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp <- vlp[vlp$uniqueDiseases >= 1, ]

vlp_summary <- tmp %>%
  mutate(
    maxMAF_bin = cut(
      maxMAF,
      breaks = maf_breaks,
      include.lowest = TRUE,
      labels = maf_labels
    )
  ) %>%
  group_by(maxMAF_bin) %>%
  summarise(
    n_total = n(),
    n_unique1 = sum(uniqueDiseases > 1, na.rm = TRUE),
    proportion_unique1 = n_unique1 / n_total,
    # compute standard error and CI for proportion
    se = sqrt(proportion_unique1 * (1 - proportion_unique1) / n_total),
    lower = proportion_unique1 - 1.96 * se,
    upper = proportion_unique1 + 1.96 * se,
    .groups = "drop"
  )

# Plot
ggplot(vlp_summary, aes(x = maxMAF_bin, y = proportion_unique1)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.002) +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Proportion of pleiotopic varaints",
    title = "Proportion of uniqueDiseases > 1 by maxMAF bin"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Confouning in each MAF bin



library(dplyr)
library(ggplot2)

# Define MAF bins and labels
maf_breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp <- vlp[vlp$uniqueDiseases >= 1, ] %>%
  mutate(
    maxMAF_bin = cut(
      maxMAF,
      breaks = maf_breaks,
      include.lowest = TRUE,
      labels = maf_labels
    )
  )

# Fit global linear model
fit_global <- lm(uniqueDiseases ~ maxCoefficientDetermination*maxEffectiveSampleSize, data = tmp)

# Add predictions
tmp <- tmp %>%
  mutate(pred = predict(fit_global, newdata = tmp))

# Compute bin-wise R2 using correlation between pred and observed
vlp_r2 <- tmp %>%
  group_by(maxMAF_bin) %>%
  summarise(
    R2 = cor(pred, uniqueDiseases, use = "complete.obs")^2,
    n = n(),
    .groups = "drop"
  )

# Plot R2 by bin
ggplot(vlp_r2, aes(x = maxMAF_bin, y = R2)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = expression(R^2),
    title = expression("Bin-wise R"^2*" from global model: uniqueDiseases ~ maxEffectiveSampleSize")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




library(dplyr)
library(ggplot2)

# Define MAF bins and labels
maf_breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

vlp1=vlp[vlp$uniqueDiseases>1,]
# Bin the data
tmp <- vlp1 %>%
  mutate(
    maxMAF_bin = cut(
      maxMAF,
      breaks = maf_breaks,
      include.lowest = TRUE,
      labels = maf_labels
    ),
    product = maxVarG*maxAbsBeta^2*maxEffectiveSampleSize
  )

# Summarise mean per bin
maf_summary <- tmp %>%
  group_by(maxMAF_bin) %>%
  summarise(
    mean_product = mean(product, na.rm = TRUE),
    se_product = sd(product, na.rm = TRUE) / sqrt(sum(!is.na(product))),
    lower = mean_product - 1.96 * se_product,
    upper = mean_product + 1.96 * se_product,
    .groups = "drop"
  )

# Plot mean ± 95% CI per bin
ggplot(maf_summary, aes(x = maxMAF_bin, y = mean_product)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.002, color = "steelblue") +
  geom_line(aes(group = 1), color = "steelblue", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "maxMAF bin",
    y = "Mean of ( maxVarG*maxAbsBeta^2*maxEffectiveSampleSize) ± 95% CI",
    title = "Mean product of  maxVarG*maxAbsBeta^2*maxEffectiveSampleSize by MAF bin for PS>1 ONLY"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





library(dplyr)
library(ggplot2)

# Prepare data
tmp <- vlp %>%
  mutate(
    # Group uniqueDiseases > 10 into "10+"
    uniqueDiseases_cat = ifelse(uniqueDiseases > 9, "9+", as.character(uniqueDiseases)),
    # Compute product of interest
    product =  maxEffectiveSampleSize*(maxAbsBeta^2)*maxVarG
  )

# Summarize mean and CI per uniqueDiseases category
summary_df <- tmp %>%
  group_by(uniqueDiseases_cat) %>%
  summarise(
    mean_product = mean(product, na.rm = TRUE),
    se_product = sd(product, na.rm = TRUE) / sqrt(sum(!is.na(product))),
    lower = mean_product - 1.96 * se_product,
    upper = mean_product + 1.96 * se_product,
    n = n(),
    .groups = "drop"
  )

# Ensure proper ordering of categories (1–9, then 9+)
summary_df$uniqueDiseases_cat <- factor(
  summary_df$uniqueDiseases_cat,
  levels = c(as.character(1:9), "9+")
)

# Plot mean ± 95% CI
ggplot(summary_df, aes(x = uniqueDiseases_cat, y = mean_product)) +
  geom_point(size = 3, color = "darkorange") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange") +
  geom_line(aes(group = 1), color = "darkorange", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Number of Unique Diseases",
    y = expression("Mean(" ~ maxVarG ~ "×" ~ maxEffectiveSampleSize ~ "×" ~ maxAbsBeta^2 ~ ") ± 95% CI"),
    title = "Mean product of maxVarG × maxEffectiveSampleSize × maxAbsBeta² by uniqueDiseases"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###########

library(dplyr)
library(ggplot2)

ind_pleiotropy=which(vlp$uniqueDiseases>1) 


# Define the bin breaks
maf_breaks <- c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)  # 0, 0.01, 0.02, …, 0.5
maf_labels <- paste0(head(maf_breaks, -1), "-", tail(maf_breaks, -1))

tmp=vlp[vlp$uniqueDiseases>=1,]

#tmp$res1=lm(uniqueDiseases~maxCoefficientDetermination*maxEffectiveSampleSize,data=tmp)$residuals
tmp$res1=tmp$res1-min(tmp$res1)+1

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
    y = "Mean PS ± 95% CI",
    title = "Mean PS by maxMAF bin for ONLY PS>1"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
