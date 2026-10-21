# Link Misspecification Study: Model M1 (Reference - Logistic Link)
# Dataset: Arrhythmia full dataset (n = 452, d = 258)
# Prior:   Regularised horseshoe prior (Piironen & Vehtari, 2017)
# Link:    Logistic (Bernoulli-logistic likelihood)


rm(list = ls())
invisible(gc())
set.seed(1904)

library(glmnet)
library(doParallel)
library(mcmcse)
library(Rcpp)

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

# Register parallel backend with up to 50 cores for glmnet
ncores <- min(10, parallel::detectCores())
registerDoParallel(cores = ncores)

source("elppd_utils.R")
Rcpp::sourceCpp("gibbs_loocv_ess.cpp")

# ---- Parameters --------------------------------------------------------------

n_sweeps   <- 1L   # ESS sweeps over theta_0..theta_D per Gibbs iteration

S_star   <- 5e3
S_max    <- 2e5
S_step   <- 5e3
S_vals   <- seq(S_star, S_max, by = S_step)

# ---- Hyper-parameters -----------------------------

p0         <- 10    # prior guess of relevant predictors (fixed, shared)
slab_scale <- 2     # s:  slab scale
slab_df    <- 4     # nu: slab degrees of freedom
interweave_tau <- TRUE  # ASIS (Yu & Meng, 2011) non-centred tau update

# ---- Data & Model Specs ------------------------------------------------------

load("data.RData")

if (file.exists("model_specs.RData")) {
  load("model_specs.RData")
}

if (exists("data_M1")) {
  model_data <- data_M1
  y <- 1 - data$y
  X <- if (exists("X_M1")) X_M1 else model_data$X
} else {
  y <- 1 - data$y
  X <- data$X
}

n <- nrow(X)
d <- ncol(X)
D <- d - 1          # number of shrunk coefficients (D = 257)
kappa <- 1:n

cat(sprintf("Setting up Model M1 (Logistic link, Full model): n = %d, d = %d\n",
            n, d, D))
cat(sprintf("Homeostasis cases (Y*=1, Normal Rhythm): %d (%.2f%%) \n",
            sum(y), 100 * mean(y)))

# ---- Global-scale hyper-parameter ------------------------------------------
mu_hat      <- mean(y)
sigma_tilde <- 1 / sqrt(mu_hat * (1 - mu_hat))
tau0        <- p0 / (D - p0) * sigma_tilde / sqrt(n)

cat(sprintf("  p0 = %.1f, sigma_tilde = %.3f, tau0 = %.5f, slab: scale = %.1f, df = %.1f\n",
            p0, sigma_tilde, tau0, slab_scale, slab_df))

# Lasso initialization via cv.glmnet with nfolds = n, parallel = TRUE
cat(sprintf("Computing Lasso initialization (nfolds = %d, parallel = TRUE)...\n", n))
cv_fit_init <- glmnet::cv.glmnet(
  x = X[, -1, drop = FALSE],
  y = y,
  family = "binomial",
  type.measure = "deviance",
  alpha = 1,
  intercept = TRUE,
  standardize = TRUE,
  nfolds = 1e2,
  parallel = TRUE
)
beta_init <- as.numeric(coef(cv_fit_init, s = "lambda.min"))

# ---- Gibbs Sampling (component-wise ESS, no warm-up) ------------------------

cat(sprintf("Running Main Gibbs LOOCV sampling for %d iterations...\n", S_max))

t_start <- proc.time()
fit_main <- gibbs_loocv_rhs_ess_cpp(
  X          = X,
  y          = y,
  theta_init = beta_init,
  tau0       = tau0,
  slab_scale = slab_scale,
  slab_df    = slab_df,
  niters     = S_max,
  link       = "logit",
  n_sweeps   = n_sweeps,
  interweave_tau = interweave_tau
)
t_elapsed    <- (proc.time() - t_start)[3]
evals_mix    <- fit_main$evals_per_coord
beta_samples <- fit_main$samples
cat(sprintf("Main Gibbs sampling complete in %.2f s. Likelihood evaluations per coordinate: %.2f\n", 
            t_elapsed, evals_mix))
cat(sprintf("Distinct folds visited: %d of %d\n", length(unique(fit_main$folds)), n))
cat(sprintf("tau: median %.4f (prior scale tau0 = %.4f); c^2: median %.3f\n",
            median(fit_main$tau), tau0, median(fit_main$c2)))

# ---- Vectorised Log-Likelihood Matrix (Logistic Link) ------------------------

compute_log_lik_mat <- function(samples, X, y) {
  eta <- tcrossprod(samples, X)
  y_mat <- matrix(y, nrow(samples), ncol(eta), byrow = TRUE)
  softplus_mat <- log1p(exp(-abs(eta))) + pmax(eta, 0)
  y_mat * eta - softplus_mat
}

cat("Computing log-likelihood matrix (Logistic link)...\n")
log_lik_mix <- compute_log_lik_mat(beta_samples, X, y)

# ---- Sequential ELPPD History -----------------------------------------------

cat("Evaluating sequential ELPPD and MCSE checkpoints...\n")
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
    tup$numerator,
    tup$denominator,
    elppd_true = NULL,
    samples    = cbind(beta_samples[1:S_curr, , drop = FALSE],
                       tau = fit_main$tau[1:S_curr],
                       c2  = fit_main$c2[1:S_curr])
  )
  history$S         <- c(history$S,         S_curr)
  history$ess_elppd <- c(history$ess_elppd, out$ess_elppd)
  history$minESS_xi <- c(history$minESS_xi, out$minESS_xi)
  history$elppd_hat <- c(history$elppd_hat, out$elppd)
  history$bias      <- c(history$bias,      out$bias)
  history$mse       <- c(history$mse,       out$mse)
  history$asymp_se  <- c(history$asymp_se,  out$asymp_se)
}

print(data.frame(
  S         = history$S,
  elppd_hat = round(history$elppd_hat, 3),
  asymp_se  = round(history$asymp_se, 4),
  ess_elppd = round(history$ess_elppd, 1),
  minESS_xi = round(history$minESS_xi, 1)
))

# ---- Save --------------------------------------------------------------------

link <- "logit"

save(history,
     evals_mix,
     t_elapsed,
     p0,
     tau0,
     slab_scale,
     slab_df,
     interweave_tau,
     n,
     d,
     link,
     file = "M1_outputs.RData")
cat("Results saved to M1_outputs.RData\n")
