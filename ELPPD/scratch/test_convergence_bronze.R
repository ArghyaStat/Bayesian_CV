# Convergence test for Bronze
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
fold_sizes <- lengths(folds)
alpha <- 1 - mean(fold_sizes) / n

log_like <- function(theta, X, y){
  mu <- X %*% theta
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

# Bronze proposal. Let's compute its mean and covariance matrix.
# log q_b = log prior + alpha * sum(log_like)
# This is a Gaussian distribution with:
V_prop <- solve(alpha * crossprod(X) / sigma2 + diag(1 / tau2, d))
m_prop <- alpha * V_prop %*% crossprod(X, y) / sigma2

for(S_iid in c(1000, 5000, 20000, 100000)){
  elppd_ests <- numeric(reps)
  var_ests <- numeric(reps)
  
  for(r in 1:reps){
    samples_b <- mvrnorm(S_iid, mu = m_prop, Sigma = V_prop)
    log_lik_b <- t(sapply(1:S_iid, function(s) log_like(samples_b[s, ], X, y)))
    
    tup_b_is <- tuple_function(
      log_lik_mat = log_lik_b,
      kappa = kappa,
      method = "Bronze",
      weight_type = "is",
      alpha = alpha
    )
    
    num <- as.matrix(tup_b_is$numerator)
    den <- as.matrix(tup_b_is$denominator)
    
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
