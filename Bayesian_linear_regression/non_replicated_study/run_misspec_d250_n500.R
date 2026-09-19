# Non-replicated study -- Misspecified model (M2): d = 250, n = 500, eta = 0.1
# Sampler: HMC-within-Gibbs on the augmented target q_mix(theta, kappa)

rm(list = ls())
gc()
set.seed(9071)

library(MASS)
library(mcmcse)
library(Rcpp)

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

source("elppd_utils.R")
source("sim_data.R")
Rcpp::sourceCpp("gibbs_loocv.cpp")

# ---- Parameters --------------------------------------------------------------

d      <- 250
n      <- 500
eta    <- 0.1
sigma2 <- 1
tau2   <- 100

S_star <- 5e3
S_max  <- 2e5
S_step <- 5e3
S_vals <- seq(S_star, S_max, by = S_step)

kappa   <- 1:n
eps_mix <- 0.48   # HMC step size for theta | fold
L_steps <- 10
n_hmc   <- 1L     # HMC transitions per Gibbs iteration

# Vectorised S x n log-likelihood matrix
compute_log_lik_mat <- function(samples, X, y, sigma2) {
  mu    <- samples %*% t(X)
  y_mat <- matrix(y, nrow = nrow(samples), ncol = ncol(mu), byrow = TRUE)
  dnorm(y_mat, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

# ---- Data --------------------------------------------------------------------

if (!file.exists("sim_data.RData")) stop("sim_data.RData not found.")
load("sim_data.RData")

y            <- sim_dat$y
dropped_col  <- sim_dat$dropped_col
R2_dropped   <- sim_dat$R2_dropped
tau2_dropped <- sim_dat$tau2_dropped
X_curr       <- sim_dat$X[, -dropped_col, drop = FALSE]

elpd_M1       <- loo_elpd(X = sim_dat$X, y = y, sigma2 = sigma2, tau2 = tau2)
elpd_M2       <- loo_elpd(X = X_curr, y = y, sigma2 = sigma2, tau2 = tau2)
elpd_diff_val <- elpd_M1 - elpd_M2
elppd_true    <- elpd_M2

cat(sprintf("LOO-ELPD  M1: %.4f  M2: %.4f  diff: %.4f\n",
            elpd_M1, elpd_M2, elpd_diff_val))

fit_lm    <- lm(y ~ X_curr - 1)
beta_init <- as.vector(coef(fit_lm))
beta_init[is.na(beta_init)] <- 0

# ---- HMC-within-Gibbs sampling -----------------------------------------------
# Each iteration: (1) draw fold kappa | theta with the Gumbel-max trick,
#                 (2) draw theta | kappa by HMC targeting p(theta | y_{-kappa}).

system.time(fit <- gibbs_loocv_cpp(X_curr, y, beta_init, sigma2, tau2,
                                   eps_mix, L_steps, S_max, n_hmc))
acc_mix     <- fit$acc
hmc_samples <- fit$samples
cat(sprintf("HMC (within Gibbs) acceptance rate: %.4f\n", acc_mix))
cat(sprintf("Distinct folds visited: %d of %d\n", length(unique(fit$folds)), n))

log_lik_mix <- compute_log_lik_mat(hmc_samples, X_curr, y, sigma2)
rm(fit)
gc()

# ---- ELPPD history -----------------------------------------------------------

history <- list(
  S         = numeric(),
  ess_elppd = numeric(),
  minESS_xi = numeric(),
  elppd_hat = numeric(),
  bias      = numeric(),
  mse       = numeric(),
  asymp_se  = numeric()
)

for (S_curr in S_vals) {
  tup <- tuple_function(log_lik_mix, kappa, S_curr)
  out <- elppd.details(
    tup$numerator, tup$denominator,
    elppd_true = elppd_true,
    samples    = hmc_samples[1:S_curr, , drop = FALSE]
  )
  history$S         <- c(history$S,         S_curr)
  history$ess_elppd <- c(history$ess_elppd, out$ess_elppd)
  history$minESS_xi <- c(history$minESS_xi, out$minESS_xi)
  history$elppd_hat <- c(history$elppd_hat, out$elppd)
  history$bias      <- c(history$bias,      out$bias)
  history$mse       <- c(history$mse,       out$mse)
  history$asymp_se  <- c(history$asymp_se,  out$asymp_se)
}

rm(log_lik_mix, hmc_samples, out, tup)
gc()

# ---- Save --------------------------------------------------------------------

save(history, acc_mix, elpd_M1, elpd_M2, elpd_diff_val,
     R2_dropped, tau2_dropped, dropped_col, d, n, eta,
     file = "M2_outputs.RData")
cat("Saved M2_outputs.RData\n")
