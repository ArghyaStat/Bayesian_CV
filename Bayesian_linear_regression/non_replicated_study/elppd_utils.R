# ELPPD utilities for Mixture IS / LOOCV.
# log_sum_exp lives in gibbs_loocv.cpp (exported via Rcpp::export).

# Mixture IS LOO weights: returns numerator and denominator matrices (S x n)
tuple_function <- function(log_lik_mat, kappa, S_curr = NULL) {

  S <- if (is.null(S_curr)) nrow(log_lik_mat) else S_curr
  n <- ncol(log_lik_mat)

  if (length(kappa) != n) stop("length(kappa) must equal ncol(log_lik_mat)")

  K           <- max(kappa)
  fold_loglik <- matrix(0, S, K)

  for (k in seq_len(K)) {
    idx              <- which(kappa == k)
    fold_loglik[, k] <- rowSums(log_lik_mat[1:S, idx, drop = FALSE])
  }

  log_inv_fold <- -fold_loglik
  log_norm     <- apply(log_inv_fold, 1, log_sum_exp)
  w_fold       <- exp(log_inv_fold - log_norm)

  w_mat <- matrix(0, S, n)
  for (k in seq_len(K)) {
    idx          <- which(kappa == k)
    w_mat[, idx] <- matrix(w_fold[, k], nrow = S, ncol = length(idx))
  }

  list(
    numerator   = exp(log_lik_mat[1:S, , drop = FALSE]) * w_mat,
    denominator = w_mat
  )
}

# ELPPD point estimate + Monte Carlo diagnostics.
# samples: S x d HMC draws, used to compute minESS_xi over theta components.
elppd.details <- function(num, den, elppd_true = NULL, samples = NULL, ...) {

  num <- as.matrix(num)
  den <- as.matrix(den)
  S   <- nrow(num)

  if (S < 2) stop("Need at least two MCMC draws.")

  mu_num    <- colMeans(num)
  mu_den    <- colMeans(den)
  mu        <- mu_num / mu_den
  elppd_hat <- if (any(mu <= 0)) -Inf else sum(log(mu))

  # Linearised scalar for MC error estimation (predictive space)
  xi <- as.vector(num %*% (1 / mu_num) + den %*% (-1 / mu_den))

  asymp_se  <- mcmcse::mcse(xi, size = NULL, r = 3,
                             method = "bm", warn = FALSE)$se^2
  post_se   <- mcmcse::mcse(xi, size = 1,    r = 1,
                             method = "bm", warn = FALSE)$se^2
  ess_elppd <- if (!is.na(asymp_se) && asymp_se > 0) {
    S * post_se / asymp_se
  } else {
    NA_real_
  }

  # Min ESS over parameter components (samples: theta / beta)
  minESS_beta <- NA_real_
  worst_idx   <- NA_integer_
  if (!is.null(samples)) {
    ess_vec     <- vapply(seq_len(ncol(samples)),
                          function(j) mcmcse::ess(samples[, j]),
                          numeric(1))
    worst_idx   <- which.min(ess_vec)
    minESS_beta <- ess_vec[worst_idx]
  }

  bias <- if (!is.null(elppd_true)) elppd_hat - elppd_true else NA_real_
  mse  <- if (!is.null(elppd_true)) bias^2 + asymp_se      else NA_real_

  list(
    elppd       = elppd_hat,
    asymp_se    = asymp_se,
    ess_elppd   = ess_elppd,
    minESS_beta = minESS_beta,
    minESS_xi   = minESS_beta,  # retained for backwards compatibility
    worst_idx   = worst_idx,
    bias        = bias,
    mse         = mse,
    xi          = xi
  )
}
