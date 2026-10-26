# Shared set-up, sampler call and sequential ELJPD evaluation

# Max-min ordering, Vecchia neighbour sets, priors and initial values
load_inputs <- function(data_file = "sim_data.RData", m = 20) {
  if (!file.exists(data_file)) stop("Run sim_data.R first.")
  load(data_file)

  ord      <- max_min(sim_dat$locs)
  locs.ord <- sim_dat$locs[ord, , drop = FALSE]
  X.ord    <- sim_dat$X[ord, , drop = FALSE]
  y.ord    <- sim_dat$y[ord]

  NNarray <- neighbor_matrix(locs.ord, m)
  storage.mode(NNarray) <- "integer"

  fit_ls <- lm.fit(X.ord, y.ord)
  b_phi  <- max(fields::rdist(locs.ord)) / 3

  list(
    y = y.ord, X = X.ord, locs = locs.ord, NNarray = NNarray,
    fold = as.integer(sim_dat$fold[ord]), K = sim_dat$K, m = m,
    V_prior = 100 * diag(ncol(X.ord)), a_sigma = 0.01, b_sigma = 0.01,
    b_phi = b_phi,
    beta_init = unname(fit_ls$coefficients),
    s2_init = mean(fit_ls$residuals^2),
    phi_init = b_phi / 4,
    r_init = 0.5,
    truth = sim_dat$truth
  )
}

run_gibbs <- function(inp, nu, niters, tuning_phi, tuning_r, fixed_fold = 0L) {
  gibbs_kfold(
    y = inp$y, X = inp$X, locs = inp$locs, NNarray = inp$NNarray,
    fold = inp$fold, K = inp$K, nu = nu,
    phi_init = inp$phi_init, r_init = inp$r_init, beta_init = inp$beta_init,
    s2_init = inp$s2_init,
    V_prior = inp$V_prior, a_sigma = inp$a_sigma, b_sigma = inp$b_sigma,
    b_phi = inp$b_phi, tuning_phi = tuning_phi, tuning_r = tuning_r,
    niters = niters, fixed_fold = fixed_fold
  )
}

# ELJPD estimate, MCSE and ESS at each checkpoint S in S_vals
eljpd_history <- function(fit, S_vals) {
  K <- ncol(fit$fold_loglik)
  history <- list(
    S = S_vals,
    eljpd_hat = numeric(length(S_vals)),
    asymp_se  = numeric(length(S_vals)),
    ess_elpd  = numeric(length(S_vals)),
    minESS_xi = numeric(length(S_vals)),
    fold_log_pred = matrix(NA_real_, length(S_vals), K,
                           dimnames = list(S_vals, paste0("fold_", seq_len(K))))
  )
  for (i in seq_along(S_vals)) {
    S   <- S_vals[i]
    tup <- tuple_function(fit$fold_loglik, S_curr = S)
    out <- eljpd.details(tup$numerator, tup$denominator, log_c = tup$log_c,
                         samples = fit$samples[seq_len(S), , drop = FALSE])
    history$eljpd_hat[i]        <- out$eljpd_hat
    history$asymp_se[i]         <- out$asymp_se
    history$ess_elpd[i]         <- out$ess_elpd
    history$minESS_xi[i]        <- out$minESS_xi
    history$fold_log_pred[i, ]  <- out$fold_log_pred
  }
  history
}
