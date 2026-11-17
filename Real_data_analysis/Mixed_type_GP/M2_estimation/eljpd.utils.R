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

# Degeneracy check for the mixture sampler and its importance weights.
# The estimator is only meaningful if the fold indicator U moves across all K folds and no
# single draw carries most of the weight in the numerator or in any fold's denominator.
# Flags: a fold never visited, or a weight vector whose largest single draw exceeds max_share.
mixture_diagnostics <- function(fold_loglik, folds = NULL, max_share = 0.10) {
  tup <- tuple_function(fold_loglik)
  S <- nrow(tup$denominator)
  K <- ncol(tup$denominator)
  kish  <- function(w) sum(w)^2 / sum(w^2)
  share <- function(w) max(w) / sum(w)

  out <- list(
    S           = S,
    fold_visits = if (is.null(folds)) rep(NA_integer_, K) else tabulate(folds, nbins = K),
    num_ess     = kish(tup$numerator[, 1]),
    num_share   = share(tup$numerator[, 1]),
    den_ess     = apply(tup$denominator, 2, kish),
    den_share   = apply(tup$denominator, 2, share)
  )
  problems <- character()
  if (!is.null(folds) && any(out$fold_visits == 0)) {
    problems <- c(problems, sprintf("fold indicator visited only %d of %d folds",
                                    sum(out$fold_visits > 0), K))
  }
  if (out$num_share > max_share) {
    problems <- c(problems, sprintf("one draw carries %.0f%% of the numerator weight",
                                    100 * out$num_share))
  }
  if (any(out$den_share > max_share)) {
    problems <- c(problems, sprintf("one draw carries more than %.0f%% of the denominator weight in %d of %d folds",
                                    100 * max_share, sum(out$den_share > max_share), K))
  }
  out$problems <- problems
  out$ok <- length(problems) == 0
  out
}

print_mixture_diagnostics <- function(d) {
  cat("---- Mixture / importance-weight diagnostics ----\n")
  cat("Fold visits            :", d$fold_visits, "\n")
  cat(sprintf("Numerator   weight ESS : %.1f of %d draws (largest single-draw share %.3f)\n",
              d$num_ess, d$S, d$num_share))
  cat("Denominator weight ESS :", round(d$den_ess, 1), "\n")
  cat("Denominator max share  :", round(d$den_share, 3), "\n")
  if (d$ok) {
    cat("No degeneracy detected.\n")
  } else {
    msg <- paste0("DEGENERATE MIXTURE ESTIMATE: ", paste(d$problems, collapse = "; "),
                  ". ELJPD and its MCSE are not reliable.")
    cat(msg, "\n")
    warning(msg, call. = FALSE, immediate. = TRUE)
  }
  invisible(d)
}
