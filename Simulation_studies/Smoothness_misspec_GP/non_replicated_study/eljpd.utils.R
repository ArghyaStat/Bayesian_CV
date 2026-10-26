# ELJPD utilities for mixture importance sampling K-fold CV

log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Mixture IS K-fold numerator and denominator matrices (S x K), stored on a shifted log scale
tuple_function <- function(fold_loglik, S_curr = NULL) {
  fl <- as.matrix(fold_loglik)

  if (!is.null(S_curr)) {
    stopifnot(
      length(S_curr) == 1L,
      is.finite(S_curr),
      S_curr >= 2L,
      S_curr <= nrow(fl),
      S_curr == as.integer(S_curr)
    )
    fl <- fl[seq_len(S_curr), , drop = FALSE]
  }

  if (nrow(fl) < 2L || !all(is.finite(fl))) {
    stop("Need at least two draws with finite fold log likelihoods.")
  }

  S <- nrow(fl)
  K <- ncol(fl)

  # Numerator and denominator in log space
  a <- apply(-fl, 1L, log_sum_exp)
  log_num <- -a
  log_den <- sweep(-fl, 1L, a, "-")

  # Shift numerator globally and each denominator column separately
  num_shift <- max(log_num)
  den_shift <- apply(log_den, 2L, max)

  num <- matrix(exp(log_num - num_shift), S, K)
  den <- exp(sweep(log_den, 2L, den_shift, "-"))

  list(
    numerator   = num,
    denominator = den,
    log_c       = num_shift - den_shift
  )
}

# ELJPD point estimate, asymptotic MCSE and ESS
eljpd.details <- function(num, den, log_c = 0, samples = NULL) {
  num <- as.matrix(num)
  den <- as.matrix(den)
  S   <- nrow(num)
  K   <- ncol(num)

  if (S < 2) stop("Need at least two MCMC draws.")

  mu_num <- colMeans(num)
  mu_den <- colMeans(den)
  fold_log_pred <- log(mu_num) - log(mu_den) + log_c
  eljpd_hat <- sum(fold_log_pred)

  # Linearised scalar for MC error estimation (delta method)
  xi <- as.vector(num %*% (1 / mu_num) - den %*% (1 / mu_den))

  asymp_se <- mcmcse::mcse(xi, size = NULL, r = 3, method = "bm", warn = FALSE)$se
  post_se  <- mcmcse::mcse(xi, size = 1,    r = 1, method = "bm", warn = FALSE)$se
  ess_elpd <- S * (post_se / asymp_se)^2

  minESS_xi <- if (!is.null(samples) && ncol(samples) > 0) {
    min(mcmcse::ess(samples, r = 3))
  } else {
    NA
  }

  list(
    eljpd_hat     = eljpd_hat,
    fold_log_pred = fold_log_pred,
    asymp_se      = asymp_se,
    ess_elpd      = ess_elpd,
    minESS_xi     = minESS_xi,
    xi            = xi
  )
}
