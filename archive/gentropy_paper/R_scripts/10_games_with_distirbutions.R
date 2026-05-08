library(MASS)

# Separate zeros and non-zeros
non_zero=vals = c(0.07883535278110519,0.030872661279000867,0.03891898167610717,
                          0.05363171690758647,0.03850553581657284,0.13520505194742938,
                          0.06502301093338424,0.07489795645306359,0.04365983747099578,
                          0.11083026574244388,0.12445843669420355,0.07684962939896704,
                          0.11580216127797223,0.22287798151967794,0.05991694265080245,
                          0.1888225775300826,0.03891898167610717,0.05491866198626455,
                          0.20340966480076852,0.2052254670613944
              )

# Fit a Gamma distribution
fit_gamma <- fitdistr(non_zero, densfun = "gamma")
fit_gamma

# Install mclust if not already
# install.packages("mclust")

library(mclust)

# Example data
X <- c(0.0788,0.0309,0.0389,0.0536,0.0385,0.1352,0.0650,0.0749,0.0437,
       0.1108,0.1245,0.0768,0.1158,0.2229,0.0599,0.1888,0.0389,0.0549,
       0.2034,0.2052)

# Fit a 2-component Gaussian mixture
gmm <- Mclust(X, G = 2)

# Check results
gmm$parameters$pro    # mixture weights
gmm$parameters$mean   # means of components
gmm$parameters$variance$sigmasq  # variances of components

# Optional: cluster assignment for each point
gmm$classification




# Load libraries
library(flexmix)
library(ggplot2)

# Your original |β| values
X <- c(0.07883535278110519,0.030872661279000867,0.03891898167610717,
       0.05363171690758647,0.03850553581657284,0.13520505194742938,
       0.06502301093338424,0.07489795645306359,0.04365983747099578,
       0.11083026574244388,0.12445843669420355,0.07684962939896704,
       0.11580216127797223,0.22287798151967794,0.05991694265080245,
       0.1888225775300826,0.03891898167610717,0.05491866198626455,
       0.20340966480076852,0.2052254670613944
)


# Compute squared values
X2 <- X^2
data <- data.frame(X2 = X2)

# Fit 2-component Gamma mixture
set.seed(123)
mix_model <- flexmix(X2 ~ 1, data = data, k = 2)

# Cluster assignments
data$cluster <- as.factor(clusters(mix_model))

# Extract Gamma parameters
params <- parameters(mix_model)
shape1 <- params["coef.(Intercept)","Comp.1"] / params["sigma","Comp.1"]^2
scale1 <- params["sigma","Comp.1"]^2
shape2 <- params["coef.(Intercept)","Comp.2"] / params["sigma","Comp.2"]^2
scale2 <- params["sigma","Comp.2"]^2


# Create density curves for plotting
x_seq <- seq(min(X2), max(X2), length.out = 500)
dens1 <- dgamma(x_seq, shape = shape1, scale = scale1)
dens2 <- dgamma(x_seq, shape = shape2, scale = scale2)
dens_df <- data.frame(x = x_seq, dens1 = dens1, dens2 = dens2)

# Plot histogram + Gamma mixture + cluster assignments
ggplot(data, aes(x = X2)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(data = dens_df, aes(x = x, y = dens1), color = "red", size = 1) +
  geom_line(data = dens_df, aes(x = x, y = dens2), color = "darkgreen", size = 1) +
  geom_point(aes(y = 0, color = cluster), position = position_jitter(height = 0.001), size = 2) +
  labs(title = "2-Component Gamma Mixture Fit on Squared |β| Values",
       x = "|β|² values",
       y = "Density",
       color = "Cluster") +
  theme_minimal()


# Posterior probabilities
post <- posterior(mix_model) # n x 2 matrix

# Mixture density over a grid
x_seq <- seq(min(X2), max(X2), length.out = 500)

dens_mix <- numeric(length(x_seq))
for (i in 1:nrow(data)) {
  # Each component density weighted by posterior for this observation
  dens_mix <- dens_mix + post[i,1]*dgamma(x_seq, shape = X2[i]/var(X2)*2, scale = var(X2)/2) +
    post[i,2]*dgamma(x_seq, shape = X2[i]/var(X2)*2, scale = var(X2)/2)
}
dens_df <- data.frame(x = x_seq, dens = dens_mix / nrow(data))

# Plot histogram + mixture
ggplot(data, aes(x = X2)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(data = dens_df, aes(x = x, y = dens), color = "red", size = 1) +
  geom_point(aes(y = 0, color = cluster), position = position_jitter(height = 0.001), size = 2) +
  labs(title = "2-Component Gamma Mixture Fit on Squared |β| Values",
       x = "|β|² values",
       y = "Density",
       color = "Cluster") +
  theme_minimal()




library(flexmix)
out=array(NA,c(1,6))

vals = c(0.07883535278110519,0.030872661279000867,0.03891898167610717,
         0.05363171690758647,0.03850553581657284,0.13520505194742938,
         0.06502301093338424,0.07489795645306359,0.04365983747099578,
         0.11083026574244388,0.12445843669420355,0.07684962939896704,
         0.11580216127797223,0.22287798151967794,0.05991694265080245,
         0.1888225775300826,0.03891898167610717,0.05491866198626455,
         0.20340966480076852,0.2052254670613944
)
vals=vals^2
X2=as.data.frame(vals)
mix1 <- flexmix(X2 ~ 1, data = data, k = 1)
mix2 <- flexmix(X2 ~ 1, data = data, k = 2)
bic1 <- BIC(mix1)
bic2 <- BIC(mix2)
aic1 <- AIC(mix1)
aic2 <- AIC(mix2)
out[1,1]=length(vals)
out[1,2]=bic1
out[1,3]=bic2
out[1,4]=aic1
out[1,5]=aic2
out[1,6]=min(sum(mix2@cluster==2)/length(vals),1-sum(mix2@cluster==2)/length(vals))


library(flexmix)
library(ggplot2)

# Your original |β| values
X <- c(0.07883535278110519,0.030872661279000867,0.03891898167610717,
       0.05363171690758647,0.03850553581657284,0.13520505194742938,
       0.06502301093338424,0.07489795645306359,0.04365983747099578,
       0.11083026574244388,0.12445843669420355,0.07684962939896704,
       0.11580216127797223,0.22287798151967794,0.05991694265080245,
       0.1888225775300826,0.03891898167610717,0.05491866198626455,
       0.20340966480076852,0.2052254670613944
)

# Compute squared values
X2 <- X^2
data <- data.frame(X2 = X2)

# Fit 2-component Gaussian mixture
set.seed(123)
mix_model <- flexmix(X2 ~ 1, data = data, k = 2, model = FLXMRglm(family="gaussian"))

# Cluster assignments
data$cluster <- as.factor(clusters(mix_model))

# Extract Gaussian parameters (means and SDs)
params <- parameters(mix_model)
mu1 <- params["coef.(Intercept)", "Comp.1"]
sd1 <- params["sigma", "Comp.1"]
mu2 <- params["coef.(Intercept)", "Comp.2"]
sd2 <- params["sigma", "Comp.2"]

# Create density curves for plotting
x_seq <- seq(min(X2), max(X2), length.out = 500)
dens1 <- dnorm(x_seq, mean = mu1, sd = sd1)
dens2 <- dnorm(x_seq, mean = mu2, sd = sd2)
dens_df <- data.frame(x = x_seq, dens1 = dens1, dens2 = dens2)

# Plot histogram + Gaussian mixture + cluster assignments
ggplot(data, aes(x = X2)) +
  geom_histogram(aes(y = ..density..), bins = 20,
                 fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(data = dens_df, aes(x = x, y = dens1),
            color = "red", size = 1) +
  geom_line(data = dens_df, aes(x = x, y = dens2),
            color = "darkgreen", size = 1) +
  geom_point(aes(y = 0, color = cluster),
             position = position_jitter(height = 0.001), size = 2) +
  labs(title = "2-Component Gaussian Mixture Fit on Squared |β| Values",
       x = "|β|² values",
       y = "Density",
       color = "Cluster") +
  theme_minimal()

# Posterior probabilities
post <- posterior(mix_model)

# Weighted mixture density
dens_mix <- post[,1][1] * dnorm(x_seq, mean = mu1, sd = sd1) +
  post[,2][1] * dnorm(x_seq, mean = mu2, sd = sd2)

dens_df_mix <- data.frame(x = x_seq, dens = dens_mix)

# Plot histogram + mixture density
ggplot(data, aes(x = X2)) +
  geom_histogram(aes(y = ..density..), bins = 20,
                 fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(data = dens_df_mix, aes(x = x, y = dens),
            color = "red", size = 1) +
  geom_point(aes(y = 0, color = cluster),
             position = position_jitter(height = 0.001), size = 2) +
  labs(title = "2-Component Gaussian Mixture Fit on Squared |β| Values",
       x = "|β|² values",
       y = "Density",
       color = "Cluster") +
  theme_minimal()

