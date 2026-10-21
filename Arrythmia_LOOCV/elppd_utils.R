# ==============================================================================
# elppd_utils.R
# Utilities for computing Leave-One-Out ELPPD under the mixture proposal q_mix
# and estimating asymptotic Monte Carlo Standard Errors (MCSE) and Effective
# Sample Size (ESS) via mcmcse.
# ==============================================================================

log_sum_exp <- function(x) {
  max_x <- max(x)
  if (is.infinite(max_x)) return(max_x)
  max_x + log(sum(exp(x - max_x)))
}

tuple_function <- function(log_lik_mat, kappa, S_curr = NULL) {
  S <- if (is.null(S_curr)) nrow(log_lik_mat) else S_curr
  n <- ncol(log_lik_mat)

  if (length(kappa) != n) {
    stop("length(kappa) must equal ncol(log_lik_mat)")
  }

  K <- max(kappa)
  fold_loglik <- matrix(0, S, K)

  for (k in seq_len(K)) {
    idx <- which(kappa == k)
    fold_loglik[, k] <- rowSums(log_lik_mat[1:S, idx, drop = FALSE])
  }

  log_inv_fold <- -fold_loglik
  log_norm <- apply(log_inv_fold, 1, log_sum_exp)
  w_fold <- exp(log_inv_fold - log_norm)

  w_mat <- matrix(0, S, n)
  for (k in seq_len(K)) {
    idx <- which(kappa == k)
    w_mat[, idx] <- matrix(w_fold[, k], nrow = S, ncol = length(idx))
  }

  list(
    numerator = exp(log_lik_mat[1:S, , drop = FALSE]) * w_mat,
    denominator = w_mat
  )
}

elppd.details <- function(num, den, elppd_true = NULL, samples = NULL) {
  num <- as.matrix(num)
  den <- as.matrix(den)
  S <- nrow(num)

  if (S < 2) {
    stop("Need at least two MCMC draws.")
  }

  mu_num <- colMeans(num)
  mu_den <- colMeans(den)
  mu <- mu_num / mu_den
  elppd_hat <- sum(log(mu))

  xi <- as.vector(num %*% (1 / mu_num) - den %*% (1 / mu_den))

  asymp_se <- mcmcse::mcse(xi, size = NULL, r = 3, method = "bm")$se
  post_se <- mcmcse::mcse(xi, size = 1, r = 1, method = "bm")$se
  ess_elppd <- S * (post_se / asymp_se)^2

  minESS_xi <- if (!is.null(samples)) {
    min(as.numeric(coda::effectiveSize(samples)))
  } else {
    NA_real_
  }

  bias <- if (!is.null(elppd_true)) elppd_hat - elppd_true else NA_real_
  mse  <- if (!is.null(elppd_true)) bias^2 + asymp_se^2 else NA_real_

  list(
    elppd = elppd_hat,
    asymp_se = asymp_se,
    ess_elppd = ess_elppd,
    minESS_xi = minESS_xi,
    bias = bias,
    mse = mse
  )
}
