
# Data simulation function 
sim_data <- function(n, d, sigma2){
  X <- matrix(rnorm(n*d), n, d)
  beta_true <- rnorm(d)
  y <- X %*% beta_true + rnorm(n, sd = sqrt(sigma2))
  list(X = X,
       y = as.vector(y),
       beta_true = beta_true)
}

# Exact tractable LOO-ELPPD for Gaussian linear model
exact_elpd_gaussian <- function(X, y, sigma2, tau2){
  
  n <- nrow(X)
  d <- ncol(X)
  
  log_ppds <- numeric(n)
  
  for(i in 1:n){
    
    # Leave-one-out data
    
    X_loo <- X[-i, , drop = FALSE]
    y_loo <- y[-i]
    
    x_test <- matrix(X[i, ], nrow = 1)
    
    y_test <- y[i]
    
    # Exact Gaussian posterior of beta | y_-i
    V_post <- solve(crossprod(X_loo) / sigma2 + diag(1 / tau2, d))
    
    m_post <- V_post %*% crossprod(X_loo, y_loo) / sigma2
    
    # Exact posterior predictive mean and variance
    
    pred_mean <- as.numeric(x_test %*% m_post)
    
    pred_var <- as.numeric(sigma2 + x_test %*% V_post %*% t(x_test))
    
    # Exact log ppd
    
    log_ppds[i] <- dnorm(x = y_test, 
                         mean = pred_mean, 
                         sd = sqrt(pred_var),
                         log = TRUE)
  }
  
  sum(log_ppds)
}

# --- Dataset 1: n = 100, d = 10 ---

set.seed(1709)


sigma2 <- 1
tau2 <- 1e2
dat <- sim_data(n = 1e2, d = 1e1, sigma2 = sigma2)
X <- dat$X
y <- dat$y
beta_true <- dat$beta_true

# Compute exact ELPPD 
exact_elppd <- exact_elpd_gaussian(X = X,
                                   y = y,
                                   sigma2 = sigma2,
                                   tau2 = tau2)

# Save to RData
save(X, y, beta_true, exact_elppd, file = "Gaussian_data.RData")
cat("exact_elppd:", exact_elppd, "\n")
