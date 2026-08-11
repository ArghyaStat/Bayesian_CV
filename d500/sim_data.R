exact_elppd_Kfold <- function(X, y, sigma2, tau2, kappa) {

  if (missing(X) || missing(y) || missing(kappa) || missing(sigma2) || missing(tau2)) {
    stop("All arguments (X, y, sigma2, tau2, kappa) must be explicitly provided.")
  }
  if (nrow(X) != length(y)) {
    stop("Matrix dimension mismatch: nrow(X) must equal the length of vector y.")
  }

  n <- nrow(X)
  d <- ncol(X)
  K <- max(kappa)
  log_ppds <- numeric(n)

  for (k in 1:K) {

    idx_test  <- which(kappa == k)
    idx_train <- which(kappa != k)

    if (length(idx_train) == 0) next

    X_train <- X[idx_train, , drop = FALSE]
    y_train <- y[idx_train]
    X_test  <- X[idx_test, , drop = FALSE]
    y_test  <- y[idx_test]

    n_train <- length(y_train)

    if (d >= n_train) {
      Omega_inv <- solve(diag(sigma2 / tau2, n_train) + X_train %*% t(X_train))
      m_post <- t(X_train) %*% (Omega_inv %*% y_train)
      pred_mean <- as.vector(X_test %*% m_post)
      U <- X_train %*% t(X_test)
      quadratic_term <- colSums((Omega_inv %*% U) * U)
      x_norm2 <- rowSums(X_test^2)
      pred_var <- sigma2 + tau2 * x_norm2 - tau2 * quadratic_term

    } else {
      V_post <- solve(crossprod(X_train) / sigma2 + diag(1 / tau2, d))
      m_post <- V_post %*% crossprod(X_train, y_train) / sigma2
      pred_mean <- as.vector(X_test %*% m_post)
      pred_var  <- sigma2 + rowSums((X_test %*% V_post) * X_test)
    }

    log_ppds[idx_test] <- dnorm(y_test, mean = pred_mean, sd = sqrt(pred_var), log = TRUE)
  }

  return(sum(log_ppds))
}


sim_data <- function(n, d, s, eta, sigma2) {
  
  if (eta <= 0 || eta >= 1) {
    stop("eta must be between 0 and 1 (exclusive)")
  }
  if (s > d) {
    stop("Number of active features (s) cannot exceed total features (d)")
  }
  
  # Calculate theoretical SNR and signal variance
  snr <- eta / (1 - eta)
  target_signal_var <- snr * sigma2
  
  # Generate X ~ N(0, 1/n)
  X <- matrix(rnorm(n * d, mean = 0, sd = 1), nrow = n, ncol = d)
  
  # Generate sparse coefficient vector
  beta_true <- rep(0, d)
  active_indices <- sample(1:d, s)
  
  # Initial random draws
  beta_sd <- sqrt(target_signal_var / s)
  beta_true[active_indices] <- rnorm(s, mean = 0, sd = beta_sd)
  
  # Exact scaling adjustment to guarantee exact target variance
  current_signal_var <- var(as.vector(X %*% beta_true))
  if (current_signal_var > 0) {
    beta_true <- beta_true * sqrt(target_signal_var / current_signal_var)
  }
  
  # Generate noise and response vector Y
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  y <- X %*% beta_true + epsilon
  
  # Return data along with underlying parameters
  return(list(
    X = X, 
    y = as.vector(y), 
    beta_true = beta_true, 
    active_indices = active_indices,
    true_snr = snr
  ))
}



# Data Generation Loop
set.seed(1709)

mydir = this.path::here()
setwd(mydir)

eta_vals <- c(0.1, 0.5, 0.9)

n <- 100
sigma2 <- 1
d <- 500
s <- min(d, floor(n/log(d)))

# Data frame to hold verification results
verification_results <- data.frame(
  True_Eta = numeric(), True_SNR = numeric(),
  Emp_Eta = numeric(), Emp_SNR = numeric()
)

for (eta in eta_vals) {
  sim_dat <- sim_data(d = d, s = s, eta = eta, n = n, sigma2 = sigma2)
  
  # Calculate empirical metrics using the corrected data
  emp_signal_var <- var(as.vector(sim_dat$X %*% sim_dat$beta_true))
  emp_snr <- emp_signal_var / sigma2
  emp_eta <- emp_snr / (1 + emp_snr)
  
  # Store metrics in dataframe
  verification_results <- rbind(
    verification_results,
    data.frame(
      True_Eta = eta,
      True_SNR = sim_dat$true_snr,
      Emp_Eta = round(emp_eta, 4),
      Emp_SNR = round(emp_snr, 4)
    )
  )
  
  # Save the simulation data
  save(sim_dat, file = sprintf("sim_data_d%d_eta%.1f.RData", d, eta))
}

print(verification_results)


