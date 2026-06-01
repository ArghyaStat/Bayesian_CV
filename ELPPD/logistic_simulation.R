# Logistic data simulation function
sim_data <- function(n, d) {
  X <- matrix(rnorm(n * d), n, d)
  beta_true <- rnorm(d)
  eta <- as.vector(X %*% beta_true)
  p <- plogis(eta)
  y <- rbinom(n, 1, p)
  list(X = X,
       y = y,
       beta_true = beta_true)
}

# Load samplers and proposals
source("proposals.R")
source("sampler.R")

# Bernoulli log-likelihood with logistic link
log_like <- function(theta, X, y) {
  eta <- as.vector(X %*% theta)
  loss <- pmax(eta, 0) + log1p(exp(-abs(eta)))
  y * eta - loss
}

grad_log_like <- function(theta, X, y) {
  eta <- as.vector(X %*% theta)
  p <- plogis(eta)
  X * as.vector(y - p)
}

# Gaussian prior
log_prior <- function(theta) {
  sum(dnorm(theta, mean = 0, sd = sqrt(tau2), log = TRUE))
}

grad_log_prior <- function(theta) {
  -theta / tau2
}


# --- Dataset Simulation ---
set.seed(1709)

n <- 1e2
d <- 1e1
tau2 <- 1e2

dat <- sim_data(n = n, d = d)
X <- dat$X
y <- dat$y
beta_true <- dat$beta_true

# Save simulated data
save(X, y, beta_true, file = "logistic_data.RData")