# Convergence test for Zanella
setwd("/Users/arghyamukherjee/Documents/GitHub/LOOCV/ELPPD")
source("proposals.R")
source("elppd_utils.R")

library(MASS)

n <- 50
d <- 10
K <- 5
reps <- 200

sigma2 <- 1
tau2 <- 1e2

set.seed(123)

sim_data <- function(n, d){
  X <- matrix(rnorm(n*d), n, d)
  beta_true <- rnorm(d)
  y <- X %*% beta_true + rnorm(n, sd = sqrt(sigma2))
  list(X = X, y = as.vector(y), beta_true = beta_true)
}

dat <- sim_data(n, d)
X <- dat$X
y <- dat$y
kappa <- sample(rep(seq_len(K), length.out = n))
folds <- split(seq_len(n), kappa)

log_like <- function(theta, X, y){
  mu <- X %*% theta
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

log_prior <- function(theta){
  sum(dnorm(theta, mean = 0, sd = sqrt(tau2), log = TRUE))
}

# Compute post mean and covariance for each fold k
V_k <- list()
m_k <- list()
log_p_y_loo <- numeric(K)

for(k in 1:K){
  idx <- which(kappa == k)
  X_loo <- X[-idx, , drop = FALSE]
  y_loo <- y[-idx]
  
  V_k[[k]] <- solve(crossprod(X_loo) / sigma2 + diag(1 / tau2, d))
  m_k[[k]] <- V_k[[k]] %*% crossprod(X_loo, y_loo) / sigma2
  
  # Compute log marginal likelihood using log p(y) = log p(y|theta) + log p(theta) - log p(theta|y)
  # evaluated at theta = 0
  log_p_theta_post <- -0.5 * d * log(2 * pi) - 0.5 * determinant(V_k[[k]], logarithm = TRUE)$modulus - 0.5 * t(m_k[[k]]) %*% solve(V_k[[k]]) %*% m_k[[k]]
  log_p_y_loo[k] <- sum(dnorm(y_loo, 0, sqrt(sigma2), log = TRUE)) + log_prior(rep(0, d)) - log_p_theta_post
}

# Zanella proposal component probabilities are proportional to 1 / p(y_loo)
log_mix_weights <- -log_p_y_loo
log_mix_weights <- log_mix_weights - log_sum_exp(log_mix_weights)
mix_weights <- exp(log_mix_weights)

# Function to draw independent samples from Zanella proposal
draw_zanella <- function(S_iid) {
  components <- sample(1:K, size = S_iid, replace = TRUE, prob = mix_weights)
  samples <- matrix(0, S_iid, d)
  for(s in 1:S_iid) {
    k <- components[s]
    samples[s, ] <- mvrnorm(1, mu = m_k[[k]], Sigma = V_k[[k]])
  }
  samples
}

for(S_iid in c(1000, 5000, 20000, 100000)){
  elppd_ests <- numeric(reps)
  var_ests <- numeric(reps)
  
  for(r in 1:reps){
    samples_z <- draw_zanella(S_iid)
    log_lik_z <- t(sapply(1:S_iid, function(s) log_like(samples_z[s, ], X, y)))
    
    tup_z_is <- tuple_function(
      log_lik_mat = log_lik_z,
      kappa = kappa,
      method = "Zanella",
      weight_type = "is"
    )
    
    num <- as.matrix(tup_z_is$numerator)
    den <- as.matrix(tup_z_is$denominator)
    
    mu_num <- colMeans(num)
    mu_den <- colMeans(den)
    
    v <- c(1 / mu_num, -1 / mu_den)
    Y <- cbind(num, den)
    xi <- as.numeric(Y %*% v)
    
    elppd_ests[r] <- sum(log(mu_num / mu_den))
    var_ests[r] <- var(xi) / S_iid
  }
  
  sample_var <- var(elppd_ests)
  mean_est_var <- mean(var_ests)
  cat("S =", S_iid, "| Sample Var =", sample_var, "| Mean Est Var =", mean_est_var, "| Ratio =", sample_var / mean_est_var, "\n")
}
