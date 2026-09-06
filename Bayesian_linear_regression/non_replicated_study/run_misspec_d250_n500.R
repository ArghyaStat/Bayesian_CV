# Non-replicated study: Misspecified model (M2)
# Configuration: d = 250, n = 500
# Sampler: Mixture IS via HMC

rm(list = ls())
gc()

library(MASS)
library(mcmcse)

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

source("elppd_utils.R")
source("proposals.R")
source("sampler.R")
source("sim_data.R")

# Simulation parameters
d      <- 250
n      <- 500
eta    <- 0.1
sigma2 <- 1
tau2   <- 100

S_star <- 5e3
S_max  <- 1e5
S_step <- 5e3
S_vals <- seq(S_star, S_max, by = S_step)

K     <- n
kappa <- 1:n

# HMC parameters
eps_mix <- 0.48
L_steps <- 20

# Vectorised log-likelihood matrix: S x n
compute_log_lik_mat <- function(samples, X, y, sigma2) {
  mu    <- samples %*% t(X)
  y_mat <- matrix(y, nrow = nrow(samples), ncol = ncol(mu), byrow = TRUE)
  dnorm(y_mat, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

# Observation-level log-likelihood and gradient
log_like <- function(theta, X, y) {
  mu <- X %*% theta
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

grad_log_like <- function(theta, X, y) {
  mu <- X %*% theta
  r  <- y - mu
  X * as.vector(r) / sigma2
}

log_prior <- function(theta) {
  sum(dnorm(theta, mean = 0, sd = sqrt(tau2), log = TRUE))
}

grad_log_prior <- function(theta) {
  -theta / tau2
}


data_file <- sprintf("sim_data_d%d.RData", d)
if (!file.exists(data_file)) {
  stop(sprintf("Data file %s not found. Please run sim_data.R first.", data_file))
}
load(data_file)

y            <- sim_dat$y
dropped_col  <- sim_dat$dropped_col
R2_dropped   <- sim_dat$R2_dropped
tau2_dropped <- sim_dat$tau2_dropped
X_curr       <- sim_dat$X[, -dropped_col, drop = FALSE]

# Exact analytical LOO-ELPD
elpd_M1      <- loo_elpd(X = sim_dat$X, y = y, sigma2 = sigma2, tau2 = tau2)
elpd_M2      <- loo_elpd(X = sim_dat$X[, -dropped_col, drop = FALSE],
                          y = y, sigma2 = sigma2, tau2 = tau2)
elpd_diff_val <- elpd_M1 - elpd_M2
elppd_true    <- elpd_M2

cat(sprintf("Exact LOO-ELPD M1: %.4f | M2: %.4f | Difference: %.4f\n",
            elpd_M1, elpd_M2, elpd_diff_val))

fit_lm    <- lm(y ~ X_curr - 1)
beta_init <- as.vector(coef(fit_lm))
beta_init[is.na(beta_init)] <- 0

fit_mix    <- leapfrog_hmc(
  epsilon = eps_mix, L = L_steps, niters = S_max,
  U = U_mix, grad_U = grad_mix, x_init = beta_init,
  X = X_curr, y = y, kappa = kappa, log_prior = log_prior,
  grad_log_prior = grad_log_prior,
  log_like = log_like, grad_log_like = grad_log_like
)
acc_mix    <- fit_mix$acc
cat(sprintf("HMC acceptance rate: %.4f\n", acc_mix))

log_lik_mix <- compute_log_lik_mat(fit_mix$samples, X_curr, y, sigma2)
rm(fit_mix)
gc()

history <- list(
  S = numeric(), ess = numeric(), elppd_hat = numeric(),
  bias = numeric(), mse = numeric(), asymp_se = numeric()
)

for (S_curr in S_vals) {
  tup <- tuple_function(log_lik_mix, kappa = kappa, S_curr = S_curr)
  out <- elppd.details(tup$numerator, tup$denominator,
                       S = S_curr, elppd_true = elppd_true)

  ess_val           <- if (!is.null(out$ess)) out$ess else out$ess_elpd
  history$S         <- c(history$S, S_curr)
  history$ess       <- c(history$ess, ess_val)
  history$elppd_hat <- c(history$elppd_hat, out$elppd)
  history$bias      <- c(history$bias, out$bias)
  history$mse       <- c(history$mse, out$mse)
  history$asymp_se  <- c(history$asymp_se, out$asymp_se)
}

rm(log_lik_mix)
gc()

last_idx <- length(history$S)
summary <- list(
  S_max     = S_max,
  ELPPD_hat = history$elppd_hat[last_idx],
  asymp_se  = history$asymp_se[last_idx],
  ESS       = history$ess[last_idx],
  bias      = history$bias[last_idx],
  mse       = history$mse[last_idx]
)

results <- list(summary = summary, acc_mix = acc_mix, history = history)

res_file <- sprintf("d%d_n%d_misspec_results.RData", d, n)
save(results, elpd_M1, elpd_M2, elpd_diff_val, R2_dropped, tau2_dropped,
     dropped_col, d, n, file = res_file)
cat(sprintf("Saved results to %s\n", res_file))
