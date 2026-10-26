# Configuration: d = 250, n = 500

# Exact Bayesian LOO-ELPD
loo_elpd <- function(X, y, sigma2, tau2) {
  
  n <- nrow(X)
  d <- ncol(X)
  Lambda <- diag(1 / tau2, d)
  M  <- crossprod(X) + sigma2 * Lambda
  M_inv  <- solve(M)

  beta_hat <- as.vector(M_inv %*% crossprod(X, y))
  e <- as.vector(y - X %*% beta_hat)
  h <- rowSums((X %*% M_inv) * X)

  pred_mean <- y - e / (1 - h)
  pred_var <- sigma2 / (1 - h)

  sum(dnorm(y, mean = pred_mean, sd = sqrt(pred_var), log = TRUE))
}

loo_elpd_diff <- function(X, y, sigma2, tau2, dropped_col) {
  
  elpd_m1 <- loo_elpd(X, y, sigma2, tau2)
  elpd_m2 <- loo_elpd(X[, -dropped_col, drop = FALSE], y, sigma2, tau2)
  elpd_m1 - elpd_m2
}

# Simulate linear regression data; the coefficient of dropped_col is set so the exact ELPD gap equals target_delta
sim_data <- function(n, d, s, eta, sigma2,
                     target_delta = 4, dropped_col = d,
                     tau2 = 100, beta_range = c(0, 10),
                     seed = NULL) {

  stopifnot(eta > 0, eta < 1, s >= 1, s <= d,
            dropped_col >= 1, dropped_col <= d)
  if (!is.null(seed)) set.seed(seed)

  # X ~ N(0, 1/n)
  X <- matrix(rnorm(n * d, mean = 0, sd = 1 / sqrt(n)), nrow = n, ncol = d)

  beta_true      <- rep(0, d)
  active_indices <- dropped_col
  
  if (s > 1) {
    
    candidates <- setdiff(seq_len(d), dropped_col)
    noise_indices <- sample(candidates, min(s - 1, length(candidates)))
    beta_true[noise_indices] <- rnorm(length(noise_indices), mean = 0, sd = 1)

    target_snr_others  <- eta / (1 - eta)
    target_signal_var_others <- target_snr_others * sigma2

    v_others  <- X[, noise_indices, drop = FALSE] %*% beta_true[noise_indices]
    current_var_others <- var(as.vector(v_others))
    
    if (current_var_others > 0) {
      beta_true[noise_indices] <- beta_true[noise_indices] *
        sqrt(target_signal_var_others / current_var_others)
    }
    
    active_indices <- c(dropped_col, noise_indices)
  }

  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))

  make_y <- function(b) {
    bt <- beta_true; bt[dropped_col] <- b
    as.vector(X %*% bt + epsilon)
  }
  
  f <- function(b) loo_elpd_diff(X, make_y(b), sigma2, tau2, dropped_col) - target_delta

  root  <- uniroot(f, interval = beta_range, extendInt = "yes")
  beta_drop <- root$root
  beta_true[dropped_col] <- beta_drop
  y <- make_y(beta_drop)

  tau2_dropped <- var(X[, dropped_col] * beta_drop)
  R2_dropped  <- tau2_dropped / (sigma2 + tau2_dropped)
  emp_signal_var_total  <- var(as.vector(X %*% beta_true))
  emp_snr_total <- emp_signal_var_total / sigma2
  emp_eta_total <- emp_snr_total / (1 + emp_snr_total)
  delta_exact <- loo_elpd_diff(X, y, sigma2, tau2, dropped_col)

  list(
    X = X, y = y, beta_true = beta_true,
    active_indices = sort(active_indices),
    dropped_col = dropped_col,
    R2_dropped  = R2_dropped,
    tau2_dropped = tau2_dropped,
    target_delta = target_delta,
    delta_exact = delta_exact,
    prior_var = tau2,
    true_eta_others = eta,
    emp_eta_total = emp_eta_total,
    emp_snr_total = emp_snr_total,
    true_snr = eta / (1 - eta),
    true_eta = eta
  )
}

# Run data generation when executed directly
if (sys.nframe() == 0) {
  set.seed(1709)
  if (requireNamespace("this.path", quietly = TRUE)) {
    setwd(this.path::here())
  }

  n <- 500
  d <- 250
  eta <- 0.1
  sigma2 <- 1
  tau2 <- 100
  s <- min(d, max(1, floor(n / log(d))))
  target_delta <- 4.0

  sim_dat <- sim_data(n = n, d = d, s = s, eta = eta, sigma2 = sigma2,
                      target_delta = target_delta, tau2 = tau2)

  save(sim_dat, file = sprintf("sim_data_d%d.RData", d))
  cat(sprintf("Generated sim_data_d%d.RData (delta_exact = %.4f)\n",
              d, sim_dat$delta_exact))
}
