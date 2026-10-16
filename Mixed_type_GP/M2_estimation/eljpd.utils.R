# ELPD utilities for Mixture IS / K-fold CV

log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Mixture IS K-fold weights: returns numerator and denominator matrices (S x K)
tuple_function <- function(fold_loglik, S_curr = NULL) {
  S <- if (is.null(S_curr)) nrow(fold_loglik) else S_curr
  K <- ncol(fold_loglik)
  fl <- fold_loglik[1:S, , drop = FALSE]

  log_inv  <- -fl
  log_norm <- apply(log_inv, 1, log_sum_exp)
  w        <- exp(log_inv - log_norm)

  log_c    <- max(-log_norm)
  u        <- matrix(exp(-log_norm - log_c), nrow = S, ncol = K)

  list(
    numerator   = u,
    denominator = w,
    log_c       = log_c
  )
}

# ELPD point estimate + Monte Carlo diagnostics
eljpd.details <- function(num, den, log_c = 0, samples = NULL) {
  num <- as.matrix(num)
  den <- as.matrix(den)
  S   <- nrow(num)
  K   <- ncol(num)

  if (S < 2) stop("Need at least two MCMC draws.")

  mu_num    <- colMeans(num)
  mu_den    <- colMeans(den)
  mu        <- mu_num / mu_den
  eljpd_hat <- sum(log(mu)) + K * log_c

  # Linearised scalar for MC error estimation
  xi <- as.vector(num %*% (1 / mu_num) - den %*% (1 / mu_den))

  asymp_se  <- mcmcse::mcse(xi, size = NULL, r = 3, method = "bm", warn = FALSE)$se
  post_se   <- mcmcse::mcse(xi, size = 1,    r = 1, method = "bm", warn = FALSE)$se
  ess_elpd  <- S * (post_se / asymp_se)^2

  minESS_xi <- if (!is.null(samples) && ncol(samples) > 0) {
    min(mcmcse::ess(samples, r = 3))
  } else {
    NA_real_
  }

  list(
    eljpd_hat = eljpd_hat,
    asymp_se  = asymp_se,
    ess_elpd  = ess_elpd,
    minESS_xi = minESS_xi,
    xi        = xi
  )
}
