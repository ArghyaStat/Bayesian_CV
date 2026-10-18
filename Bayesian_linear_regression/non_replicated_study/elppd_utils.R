# ELPPD utilities for Mixture Importance Sampling.

# Stable log-sum-exp
log_sum_exp <- function(x) {
  if (any(is.infinite(x) & x > 0)) return(Inf)
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}

# Mixture IS LOO weights 
tuple_function <- function(log_lik_mat, kappa, S_curr = NULL) {
  if (is.null(S_curr)) {
    S <- nrow(log_lik_mat)
  } else {
    S <- S_curr
  }
  n <- ncol(log_lik_mat)

  if (length(kappa) != n) {
    stop("length(kappa) must equal ncol(log_lik_mat)")
  }

  K <- max(kappa)

  # Fold log-likelihoods: S x K
  fold_loglik <- matrix(0, S, K)
  for (k in seq_len(K)) {
    idx <- which(kappa == k)
    fold_loglik[, k] <- rowSums(log_lik_mat[1:S, idx, drop = FALSE])
  }

  # Mixture weights
  log_inv_fold <- - fold_loglik
  log_norm <- apply(log_inv_fold, 1, log_sum_exp)
  log_w <- log_inv_fold - log_norm
  w_fold <- exp(log_w)

  # Map fold weights to observations
  w_mat <- matrix(0, S, n)
  for (k in seq_len(K)) {
    idx <- which(kappa == k)
    w_mat[, idx] <- matrix(w_fold[, k], nrow = S, ncol = length(idx))
  }

  numerator   <- exp(log_lik_mat[1:S, , drop = FALSE]) * w_mat
  denominator <- w_mat

  list(numerator = numerator, denominator = denominator)
}

# Asymptotic variance estimator and ELPPD calculation
elppd.details <- function(num, den, S = NULL, b = NULL,
                          elppd_true = NULL, ...) {
  num <- as.matrix(num)
  den <- as.matrix(den)

  mu_num <- colMeans(num)
  mu_den <- colMeans(den)

  S <- nrow(num)
  n <- ncol(num)

  if (S < 2) {
    stop("Need at least two MCMC draws.")
  }

  se2_num <- numeric(n)
  se2_den <- numeric(n)

  for (i in seq_len(n)) {
    se2_num[i] <- mcmcse::mcse(num[, i])$se^2
    se2_den[i] <- mcmcse::mcse(den[, i])$se^2
  }

  v_num <- mu_num / (mu_num^2 + se2_num)
  v_den <- - mu_den / (mu_den^2 + se2_den)

  # ELPPD point estimate
  mu        <- mu_num / mu_den
  elppd_hat <- sum(log(mu))

  if (any(mu <= 0)) {
    elppd_hat <- -Inf
  }

  # Batch means Monte Carlo error
  xi <- as.vector(num %*% v_num + den %*% v_den)

  asymp_se <- (mcmcse::mcse(xi, size = NULL, r = 3,
                            method = "bm", warn = FALSE)$se^2)
  post_se  <- (mcmcse::mcse(xi, size = 1, r = 3,
                            method = "bm", warn = FALSE)$se^2)

  ess <- if (!is.na(asymp_se) && asymp_se > 0) {
    S * (post_se / asymp_se)
  } else NA_real_

  if (!is.null(elppd_true)) {
    bias <- elppd_hat - elppd_true
    mse  <- (bias^2) + asymp_se
  } else {
    bias <- NA_real_
    mse  <- NA_real_
  }

  list(elppd    = elppd_hat,
       asymp_se = asymp_se,
       ess      = ess,
       bias     = bias,
       mse      = mse)
}
