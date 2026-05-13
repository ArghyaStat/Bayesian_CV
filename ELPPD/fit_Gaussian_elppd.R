rm(list = ls())

mydir <- this.path::here()
setwd(mydir)


library(MASS)
library(doParallel)
library(foreach)

source("proposals.R")
source("sampler.R")
source("elppd_components.R")


# Simulation design


n <- 50
p <- 10
K <- 5
S <- 1e4

sigma2 <- 1
tau2 <- 1e2

set.seed(1709)


# Simulate Gaussian linear model

sim_data <- function(n, p){
  
  X <- matrix(rnorm(n*p), n, p)
  
  beta_true <- rnorm(p)
  
  y <- X %*% beta_true + rnorm(n, sd = sqrt(sigma2))
  
  list(
    X = X,
    y = as.vector(y),
    beta_true = beta_true
  )
}

dat <- sim_data(n, p)

X <- dat$X
y <- dat$y

reps <- 100
n.cores <- 5

# Gaussian likelihood

log_like <- function(theta, X, y){
  
  mu <- X %*% theta
  
  dnorm(y, 
        mean = mu,
        sd = sqrt(sigma2),
        log = TRUE)
}


# Gradient log-likelihood

grad_log_like <- function(theta, X, y){
  
  mu <- X %*% theta
  
  r <- y - mu
  
  X * as.vector(r) / sigma2
}

# Gaussian prior

log_prior <- function(theta){
  
  sum(dnorm(theta,
            mean = 0,
            sd = sqrt(tau2),
            log = TRUE))
}


# Gradient prior

grad_log_prior <- function(theta){
  
  -theta / tau2
}

# Parallel for each loop

cl <- makeCluster(n.cores)

registerDoParallel(cl)

# Main simulation loop

results <- foreach(r = 1:reps, 
                   .combine = rbind) %dopar% {
  
 
  # Fold construction
  
  fold_id <- sample(rep(seq_len(K), length.out = n))
  
  folds <- split(seq_len(n), fold_id)
  
  fold_sizes <- lengths(folds)

  # Fold allocation map kappa(i)
  
  kappa <- integer(n)
  
  for(k in seq_len(K)){
    
    kappa[folds[[k]]] <- k
  }
  
 
  # Fold-specific powers alpha_k

  alpha_fold <- (n - fold_sizes) / n
 
  # Scalar Bronze power
  alpha <- mean(alpha_fold)
  
  # Initial value
  x_init <- lm(y ~ X - 1)$coef
  
  # Gelfand sampler
  
   fit_g <- leapfrog_hmc(
    epsilon = 0.16,
    L = 10,
    niters = S,
    U = U_g,
    grad_U = grad_g,
    x_init = x_init,
    M = diag(p),
    X = X,
    y = y,
    log_like = log_like,
    log_prior = log_prior,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior)
  

  # Bronze sampler
  
  fit_b <- leapfrog_hmc(
    epsilon = 0.16,
    L = 10,
    niters = S,
    U = U_b,
    grad_U = grad_b,
    x_init = x_init,
    M = diag(p),
    X = X,
    y = y,
    alpha = alpha,
    log_like = log_like,
    log_prior = log_prior,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior
  )
  

  # Zanella sampler
  
  fit_z <- leapfrog_hmc(
    epsilon = 0.15,
    L = 10,
    niters = S,
    U = U_z,
    grad_U = grad_z,
    x_init = x_init,
    M = diag(p),
    X = X,
    y = y,
    log_like = log_like,
    log_prior = log_prior,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior)
  
  samples_g <- fit_g$samples
  samples_b <- fit_b$samples
  samples_z <- fit_z$samples
  
  # Log-likelihood matrices of dimension S \times n
  
  log_lik_g <- t(sapply(1:S, function(s) 
    log_like(samples_g[s, ], X, y)))
  
  
  log_lik_b <- t(sapply(1:S, function(s) 
    log_like(samples_b[s, ], X, y)))
  
  
  log_lik_z <- t(sapply(1:S, function(s) 
    log_like(samples_z[s, ], X, y)))
  
  # Tuple construction
  
  tup_g <- tuple_function(log_lik_mat = log_lik_g, 
                          method = "Gelfand")
 
  tup_b <- tuple_function(log_lik_mat = log_lik_b,
                          method = "Bronze",
                          alpha = alpha)

  tup_z <- tuple_function(log_lik_mat = log_lik_z,
                          method = "Zanella")
  
  
  # Asymptotic variance estimation
  
  out_g <- asymp.var.elppd(num = tup_g$numerator,
                           den = tup_g$denominator,
                           S = S)
  
  out_b <- asymp.var.elppd(num = tup_b$numerator,
                           den = tup_b$denominator,
                           S = S)
  
  out_z <- asymp.var.elppd(num = tup_z$numerator,
                           den = tup_z$denominator,
                           S = S)
  
  
  c(elppd_g = out_g$elppd,
    var_g = out_g$asymp_var,
    acc_g = fit_g$acc,
    
    elppd_b = out_b$elppd,
    var_b = out_b$asymp_var,
    acc_b = fit_b$acc,
    
    elppd_z = out_z$elppd,
    var_z = out_z$asymp_var,
    acc_z = fit_z$acc
  )
}

stopCluster(cl)

results <- as.data.frame(results)

save(results,file = paste0(K,
    "_fold_elppd_results.RData"))