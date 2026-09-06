rm(list = ls())
gc()

library(MASS)
library(foreach)
library(doParallel)
library(loo)
library(mcmcse)
library(glmnet)

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

source("elppd_utils.R")
source("proposals.R")
source("sampler.R")
source("sim_data.R")

############################################################
######### Simulation Parameters ############################
############################################################
d <- 500
eta <- 0.1
sigma2 <- 1
tau2 <- 100/d
n <- 100
reps <- 100

S_star <- 4e4
S_max <- 1e5
S_step <- 5e3
nom_width <- 0.05
z_crit <- qnorm(0.975)
ess_target <- 4 * z_crit^2 / (nom_width^2)
K.vec <- c(5, 10, 20, 50, 100)

# Load simulated data
load(sprintf("sim_data_d%d_eta%.1f.RData", d, eta))
X <- sim_dat$X
y <- sim_dat$y
beta_true <- sim_dat$beta_true

log_like <- function(theta, X, y) {
  mu <- X %*% theta
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

grad_log_like <- function(theta, X, y) {
  mu <- X %*% theta
  r <- y - mu
  X * as.vector(r) / sigma2
}

log_prior <- function(theta) {
  sum(dnorm(theta, mean = 0, sd = sqrt(tau2), log = TRUE))
}

grad_log_prior <- function(theta) {
  -theta / tau2
}

# Model-specific helper functions for potential and gradient
U_fold <- function(theta, X_tr, y_tr) {
  lp <- log_prior(theta)
  mu_tr <- as.vector(X_tr %*% theta)
  ll <- sum(dnorm(y_tr, mean = mu_tr, sd = sqrt(sigma2), log = TRUE))
  return(-(lp + ll))
}

grad_U_fold <- function(theta, X_tr, y_tr) {
  grad_lp <- grad_log_prior(theta)
  r <- y_tr - as.vector(X_tr %*% theta)
  grad_ll <- as.vector(crossprod(X_tr, r)) / sigma2
  return(-(grad_lp + grad_ll))
}

# Initialize parameter starting values
if (d < n) {
  fit_lm <- lm(y ~ X - 1)
  beta_init <- as.vector(coef(fit_lm))
  beta_init[is.na(beta_init)] <- 0
} else {
  fit_lasso <- glmnet::cv.glmnet(X, y, standardise = FALSE, intercept = FALSE)
  beta_init <- as.vector(coef(fit_lasso, s = "lambda.min"))[-1]
}

set.seed(2026)
eps_post <- 7e-2
eps_b <- 7e-2
eps_mix <- 3e-2
K = 5

fold_id <- sample(rep(seq_len(K), length.out = n))
folds <- split(seq_len(n), fold_id)
fold_sizes <- lengths(folds)

kappa <- integer(n)
for (k in seq_len(K)) {
  kappa[folds[[k]]] <- k
}

elppd_lin_reg <- exact_elppd_Kfold(X = X, y = y, 
                                   sigma2 = sigma2, 
                                   tau2 = tau2, 
                                   kappa = kappa)

alpha_fold <- (n - fold_sizes) / n
alpha <- mean(alpha_fold)

# Setup history storage
est_keys <- c("post_is", "post_tis", "post_psis",
              "b_is", "b_tis", "b_psis",
              "mix_is", "mix_tis", "mix_psis")

sterm_res <- rep(NA_real_, 9)
names(sterm_res) <- est_keys

ess_res <- rep(NA_real_, 9)
names(ess_res) <- est_keys

history <- list()
for (key in est_keys) {
  history[[key]] <- list(
    S = numeric(0),
    ess = numeric(0),
    k_exceeding_0.7 = list(),
    elppd_hat = numeric(0),
    bias = numeric(0),
    mse = numeric(0),
    S_term = NA_integer_,
    asymp_se = numeric(0)
  )
}

S_vals <- seq(S_star, S_max, by = S_step)

# --- Method 1: Posterior proposal ---
fit_post_start <- Sys.time()
fit_post <- leapfrog_hmc(
  epsilon = eps_post, L = 10, niters = S_max,
  U = U_post, grad_U = grad_post, x_init = beta_init,
  X = X, y = y, log_prior = log_prior, grad_log_prior = grad_log_prior,
  log_like = log_like, grad_log_like = grad_log_like
)
acc_post <- fit_post$acc
samples_post <- fit_post$samples
beta.mean.post <- colMeans(samples_post)
rm(fit_post)
gc()

log_lik_post <- matrix(0, S_max, n)
for (s in 1:S_max) {
  log_lik_post[s, ] <- log_like(samples_post[s, ], X, y)
}

for (S_curr in S_vals) {
  for (w_idx in 1:3) {
    e <- w_idx
    w_type <- c("is", "tis", "psis")[w_idx]
    tup <- tuple_function(log_lik_post, kappa, method = "posterior",
                          weight_type = w_type, alpha = alpha,
                          S_curr = S_curr)
    out <- elppd.details(tup$numerator, tup$denominator, S = S_curr,
                         b = NULL, method = "mcmcse", elppd_true = elppd_lin_reg)
    history[[e]]$S <- c(history[[e]]$S, S_curr)
    history[[e]]$ess <- c(history[[e]]$ess, out$ess)
    history[[e]]$k_exceeding_0.7[[length(history[[e]]$k_exceeding_0.7) + 1]] <- tup$pct_k
    history[[e]]$elppd_hat <- c(history[[e]]$elppd_hat, out$elppd)
    history[[e]]$bias <- c(history[[e]]$bias, out$bias)
    history[[e]]$mse <- c(history[[e]]$mse, out$mse)
    history[[e]]$asymp_se <- c(history[[e]]$asymp_se, out$asymp_se)
    
    if (is.na(sterm_res[e]) && S_curr >= S_star &&
        !is.na(out$ess) && out$ess >= ess_target) {
      sterm_res[e] <- S_curr
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_curr
    }
    if (S_curr == S_max && is.na(sterm_res[e])) {
      sterm_res[e] <- S_max
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_max
    }
  }
}
rm(log_lik_post)
gc()
fit_post_end <- Sys.time()
fit_post_time <- as.numeric(difftime(fit_post_end, fit_post_start, units = "mins"))

# Cheap MAP and SVD coloring variables setup using HMC samples
modes <- matrix(NA_real_, nrow = K, ncol = d)
P_list <- list()
B_list <- list()
B_prime_list <- list()
log_det_V <- numeric(K)

m <- colMeans(samples_post)
idx_sub <- seq(1, S_max, length.out = min(S_max, 2000))
V <- cov(samples_post[idx_sub, ]) + diag(1e-6, d)
d_vec <- pmax(diag(V), 1e-6)

grad_m <- grad_log_like(m, X, y)

for (k in seq_len(K)) {
  idx_test <- which(kappa == k)
  G_Sk <- grad_m[idx_test, , drop = FALSE]
  temp <- G_Sk %*% V
  inner_inv <- solve(diag(length(idx_test)) - temp %*% t(G_Sk))
  V_k <- V + t(temp) %*% inner_inv %*% temp
  
  modes[k, ] <- m + as.vector(V_k %*% colSums(G_Sk))
  
  C_k <- diag(1, d) - (V_k / sqrt(d_vec)) / rep(sqrt(d_vec), each = d)
  eig_C <- eigen(C_k)
  r <- length(idx_test)
  Sigma_sq <- eig_C$values[1:r]
  Sigma <- sqrt(pmax(pmin(Sigma_sq, 0.9999), 0))
  
  P_list[[k]] <- eig_C$vectors[, 1:r, drop = FALSE]
  B_list[[k]] <- sqrt(1 - Sigma^2) - 1
  B_prime_list[[k]] <- 1 / sqrt(1 - Sigma^2) - 1
  log_det_V[k] <- sum(log(1 - Sigma^2))
}
rm(samples_post)
gc()

# --- Method 3: Mixture proposal ---
fit_mix_start <- Sys.time()
fit_mix <- non_rev_sampler(
  epsilon = eps_mix, L = 10, xi = 1, niters = S_max,
  x_init = beta_init, M = diag(d),
  X = X, y = y, log_prior = log_prior, grad_log_prior = grad_log_prior,
  log_like = log_like, grad_log_like = grad_log_like,
  kappa = kappa, modes = modes, P_list = P_list, B_list = B_list, B_prime_list = B_prime_list,
  log_det_V = log_det_V, v_prior = d_vec
)
acc_mix <- fit_mix$acc
samples_mix <- fit_mix$samples
beta.mean.mix <- colMeans(samples_mix)
rm(fit_mix)
gc()

log_lik_mix <- matrix(0, S_max, n)
for (s in 1:S_max) {
  log_lik_mix[s, ] <- log_like(samples_mix[s, ], X, y)
}
rm(samples_mix)
gc()

for (S_curr in S_vals) {
  for (w_idx in 1:3) {
    e <- 6 + w_idx
    w_type <- c("is", "tis", "psis")[w_idx]
    tup <- tuple_function(log_lik_mix, kappa, method = "mixture",
                          weight_type = w_type, alpha = alpha,
                          S_curr = S_curr)
    out <- elppd.details(tup$numerator, tup$denominator, S = S_curr,
                         b = NULL, method = "mcmcse", elppd_true = elppd_lin_reg)
    history[[e]]$S <- c(history[[e]]$S, S_curr)
    history[[e]]$ess <- c(history[[e]]$ess, out$ess)
    history[[e]]$k_exceeding_0.7[[length(history[[e]]$k_exceeding_0.7) + 1]] <- tup$pct_k
    history[[e]]$elppd_hat <- c(history[[e]]$elppd_hat, out$elppd)
    history[[e]]$bias <- c(history[[e]]$bias, out$bias)
    history[[e]]$mse <- c(history[[e]]$mse, out$mse)
    history[[e]]$asymp_se <- c(history[[e]]$asymp_se, out$asymp_se)
    
    if (is.na(sterm_res[e]) && S_curr >= S_star &&
        !is.na(out$ess) && out$ess >= ess_target) {
      sterm_res[e] <- S_curr
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_curr
    }
    if (S_curr == S_max && is.na(sterm_res[e])) {
      sterm_res[e] <- S_max
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_max
    }
  }
}
rm(log_lik_mix)
gc()
fit_mix_end <- Sys.time()
fit_mix_time <- as.numeric(difftime(fit_mix_end, fit_mix_start, units = "mins"))

# --- Method 2: Bronze proposal ---
fit_b_start <- Sys.time()
fit_b <- leapfrog_hmc(
  epsilon = eps_b, L = 10, niters = S_max,
  U = U_b, grad_U = grad_b, x_init = beta_init,
  X = X, y = y, log_prior = log_prior, grad_log_prior = grad_log_prior,
  log_like = log_like, grad_log_like = grad_log_like,
  alpha = alpha
)
acc_b <- fit_b$acc
samples_b <- fit_b$samples
beta.mean.b <- colMeans(samples_b)
rm(fit_b)
gc()

log_lik_b <- matrix(0, S_max, n)
for (s in 1:S_max) {
  log_lik_b[s, ] <- log_like(samples_b[s, ], X, y)
}
rm(samples_b)
gc()

for (S_curr in S_vals) {
  for (w_idx in 1:3) {
    e <- 3 + w_idx
    w_type <- c("is", "tis", "psis")[w_idx]
    tup <- tuple_function(log_lik_b, kappa, method = "Bronze",
                          weight_type = w_type, alpha = alpha,
                          S_curr = S_curr)
    out <- elppd.details(tup$numerator, tup$denominator, S = S_curr,
                         b = NULL, method = "mcmcse", elppd_true = elppd_lin_reg)
    history[[e]]$S <- c(history[[e]]$S, S_curr)
    history[[e]]$ess <- c(history[[e]]$ess, out$ess)
    history[[e]]$k_exceeding_0.7[[length(history[[e]]$k_exceeding_0.7) + 1]] <- tup$pct_k
    history[[e]]$elppd_hat <- c(history[[e]]$elppd_hat, out$elppd)
    history[[e]]$bias <- c(history[[e]]$bias, out$bias)
    history[[e]]$mse <- c(history[[e]]$mse, out$mse)
    history[[e]]$asymp_se <- c(history[[e]]$asymp_se, out$asymp_se)
    
    if (is.na(sterm_res[e]) && S_curr >= S_star &&
        !is.na(out$ess) && out$ess >= ess_target) {
      sterm_res[e] <- S_curr
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_curr
    }
    if (S_curr == S_max && is.na(sterm_res[e])) {
      sterm_res[e] <- S_max
      ess_res[e] <- out$ess
      history[[e]]$S_term <- S_max
    }
  }
}
rm(log_lik_b)
gc()
fit_b_end <- Sys.time()
fit_b_time <- as.numeric(difftime(fit_b_end, fit_b_start, units = "mins"))

result <- list(
  history = history,
  sterm_res = sterm_res,
  ess_res = ess_res,
  elppd_lin_reg = elppd_lin_reg,
  acc_post = acc_post,
  acc_b = acc_b,
  acc_mix = acc_mix,
  fit_mix_time = fit_mix_time,
  fit_post_time = fit_post_time,
  fit_b_time = fit_b_time,
  beta.mean.mix = beta.mean.mix,
  beta.mean.b = beta.mean.b,
  beta.mean.post = beta.mean.post
)

res_file <- sprintf("d%d_eta%.1f_K%d_fold_elppd_results_single.RData", d, eta, K)
save(result, file = res_file)
cat(sprintf("Saved: %s\n", res_file))
