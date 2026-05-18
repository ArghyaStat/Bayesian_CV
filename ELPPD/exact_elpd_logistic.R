rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

library(MASS)
library(doParallel)
library(foreach)

source("proposals.R")
source("sampler.R")
source("elppd_utils.R")


n <- 50
d <- 10
S <- 5e4

# Hyperparameter

tau2 <- 1e2

set.seed(1709)


# Simulate Logistic Linear Model

sim_data <- function(n, d){
  
  X <- matrix(rnorm(n * d), n, d)
  
  beta_true <- rnorm(d)
  
  eta <- X %*% beta_true
  
  prob <- plogis(eta)
  
  y <- rbinom(
    n,
    size = 1,
    prob = prob
  )
  
  list(
    X = X,
    y = y,
    beta_true = beta_true
  )
}

dat <- sim_data(n, d)

X <- dat$X
y <- dat$y



log_like <- function(theta, X, y){
  
  eta <- as.vector(X %*% theta)
  
  out <- numeric(length(eta))
  
  idx <- eta > 0
  
  out[idx] <-
    y[idx] * eta[idx] -
    (eta[idx] + log1p(exp(-eta[idx])))
  
  out[!idx] <-
    y[!idx] * eta[!idx] -
    log1p(exp(eta[!idx]))
  
  out
}


grad_log_like <- function(theta, X, y){
  
  eta <- as.vector(X %*% theta)
  
  prob <- plogis(eta)
  
  X * as.vector(y - prob)
}


log_prior <- function(theta){
  
  sum(
    dnorm(
      theta,
      mean = 0,
      sd = sqrt(tau2),
      log = TRUE
    )
  )
}



grad_log_prior <- function(theta){
  
  -theta / tau2
}



n.cores <- n

cl <- makeCluster(n.cores)

registerDoParallel(cl)

clusterExport(
  cl,
  c(
    "X",
    "y",
    "d",
    "S",
    "tau2",
    "log_like",
    "log_prior",
    "grad_log_like",
    "grad_log_prior",
    "leapfrog_hmc"
  ),
  envir = environment()
)


# Exact Leave-One-Out ELPPD


log_ppds <- foreach(
  i = 1:n,
  .combine = c
) %dopar% {
  
  
  # Leave-one-out data
 
  
  X_loo <- X[-i, , drop = FALSE]
  y_loo <- y[-i]
  
  x_test <- matrix(
    X[i, ],
    nrow = 1
  )
  
  y_test <- y[i]
  
 
  
  # LOO posterior target
  
  
  U_loo <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior,
                    grad_log_like = NULL,
                    grad_log_prior = NULL){
    
    -(sum(
      log_like(
        theta = theta,
        X = X_loo,
        y = y_loo
      )
    ) +
      log_prior(theta))
  }
  
  
  grad_loo <- function(theta,
                       X,
                       y,
                       log_like = NULL,
                       log_prior = NULL,
                       grad_log_like,
                       grad_log_prior){
    
    -(colSums(
      grad_log_like(
        theta = theta,
        X = X_loo,
        y = y_loo
      )
    ) +
      grad_log_prior(theta))
  }
  
 
  # Draw LOO posterior
  
  fit_loo <- leapfrog_hmc(
    epsilon = 0.56,
    L = 10,
    niters = S,
    U = U_loo,
    grad_U = grad_loo,
    x_init = rep(0, d),
    M = diag(d),
    X = X_loo,
    y = y_loo,
    log_like = log_like,
    log_prior = log_prior,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior
  )
  
 
  beta_draws <- as.matrix(fit_loo$samples)
  
  
  eta_test <- as.vector(
    beta_draws %*% t(x_test)
  )
  
  prob_test <- plogis(eta_test)
  
  log_lik_test <- dbinom(
    x = y_test,
    size = 1,
    prob = prob_test,
    log = TRUE
  )
  
  
  m <- max(log_lik_test)
  
  m + log(
    mean(
      exp(log_lik_test - m)
    )
  )
}

stopCluster(cl)


elppd_exact <- sum(log_ppds)

print(elppd_exact)

save(
  elppd_exact,
  file = "elppd_exact_logistic.RData"
)