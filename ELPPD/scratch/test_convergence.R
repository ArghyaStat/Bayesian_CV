# Convergence test for delta method
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

V_post <- solve(crossprod(X) / sigma2 + diag(1 / tau2, d))
m_post <- V_post %*% crossprod(X, y) / sigma2

for(S_iid in c(1000, 5000, 20000, 100000)){
  elppd_ests <- numeric(reps)
  var_ests <- numeric(reps)
  
  for(r in 1:reps){
    samples_g <- mvrnorm(S_iid, mu = m_post, Sigma = V_post)
    log_lik_g <- t(sapply(1:S_iid, function(s) log_like(samples_g[s, ], X, y)))
    
    tup_g_is <- tuple_function(
      log_lik_mat = log_lik_g,
      kappa = kappa,
      method = "Gelfand",
      weight_type = "is"
    )
    
    num <- as.matrix(tup_g_is$numerator)
    den <- as.matrix(tup_g_is$denominator)
    
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
