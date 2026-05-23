rm(list = ls())

mydir <- this.path::here()
setwd(mydir)


library(MASS)
library(foreach)
library(doParallel)
library(loo)

source("proposals.R")
source("sampler.R")
source("elppd_utils.R")


# Simulation design


n <- 50
d <- 10
K <- 50
S <- 1e5

# Hyperparameters

sigma2 <- 1
tau2 <- 1e2

set.seed(1709)


# Simulate Gaussian linear model

sim_data <- function(n, d){
  
  X <- matrix(rnorm(n*d), n, d)
  
  beta_true <- rnorm(d)
  
  y <- X %*% beta_true + rnorm(n, sd = sqrt(sigma2))
  
  list(X = X,
       y = as.vector(y),
       beta_true = beta_true)
}

dat <- sim_data(n, d)

X <- dat$X
y <- dat$y

# Fold construction

kappa <- sample(rep(seq_len(K),
                    length.out = n))

folds <- split(seq_len(n), kappa)

# size of each fold
fold_sizes <- lengths(folds)

# Fold-specific powers alpha_k
alpha_fold <- (n - fold_sizes) / n

alpha <- 1 - mean(fold_sizes) / n

reps <- 100
n.cores <- 50

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
                     
                     
                     # Initial value
                     theta_init <- lm(y ~ X - 1)$coef
                     
                     # Gelfand sampler
                     
                     fit_g <- leapfrog_hmc(
                       epsilon = 0.16,
                       L = 10,
                       niters = S,
                       U = U_g,
                       grad_U = grad_g,
                       x_init = theta_init,
                       M = diag(d),
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
                       x_init = theta_init,
                       M = diag(d),
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
                       epsilon = 0.16,
                       L = 10,
                       niters = S,
                       U = U_z,
                       grad_U = grad_z,
                       x_init = theta_init,
                       M = diag(d),
                       X = X,
                       y = y,
                       kappa = kappa,
                       log_like = log_like,
                       log_prior = log_prior,
                       grad_log_like = grad_log_like,
                       grad_log_prior = grad_log_prior
                     )
                     
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
                     
                     
                     # Gelfand
                     
                     tup_g_is <- tuple_function(
                       log_lik_mat = log_lik_g,
                       kappa = kappa,
                       method = "Gelfand",
                       weight_type = "is"
                     )
                     
                     tup_g_tis <- tuple_function(
                       log_lik_mat = log_lik_g,
                       kappa = kappa,
                       method = "Gelfand",
                       weight_type = "tis"
                     )
                     
                     tup_g_psis <- tuple_function(
                       log_lik_mat = log_lik_g,
                       kappa = kappa,
                       method = "Gelfand",
                       weight_type = "psis"
                     )
                     
                     
                     # Bronze
                
                     
                     tup_b_is <- tuple_function(
                       log_lik_mat = log_lik_b,
                       kappa = kappa,
                       method = "Bronze",
                       weight_type = "is",
                       alpha = alpha
                     )
                     
                     tup_b_tis <- tuple_function(
                       log_lik_mat = log_lik_b,
                       kappa = kappa,
                       method = "Bronze",
                       weight_type = "tis",
                       alpha = alpha
                     )
                     
                     tup_b_psis <- tuple_function(
                       log_lik_mat = log_lik_b,
                       kappa = kappa,
                       method = "Bronze",
                       weight_type = "psis",
                       alpha = alpha
                     )
                     
                     
                     # Zanella
                     
                     tup_z_is <- tuple_function(
                       log_lik_mat = log_lik_z,
                       kappa = kappa,
                       method = "Zanella",
                       weight_type = "is"
                     )
                     
                     tup_z_tis <- tuple_function(
                       log_lik_mat = log_lik_z,
                       kappa = kappa,
                       method = "Zanella",
                       weight_type = "tis"
                     )
                     
                     tup_z_psis <- tuple_function(
                       log_lik_mat = log_lik_z,
                       kappa = kappa,
                       method = "Zanella",
                       weight_type = "psis"
                     )
                     
                     
                     # Asymptotic SE estimation
                     
                     # Gelfand outputs
                     
                     out_g_is <- elppd.details(num = tup_g_is$numerator,
                                               den = tup_g_is$denominator,
                                               S = NULL,
                                               b = NULL,
                                               method = "mcmcse")
                     
                     
                     out_g_tis <- elppd.details(num = tup_g_tis$numerator,
                                                den = tup_g_tis$denominator,
                                                S = NULL,
                                                b = NULL,
                                                method = "mcmcse")
                     
                     
                     out_g_psis <- elppd.details(num = tup_g_psis$numerator,
                                                 den = tup_g_psis$denominator,
                                                 S = NULL,
                                                 b = NULL,
                                                 method = "mcmcse")
                     
                     # Bronze outputs
                     
                     out_b_is <- elppd.details(num = tup_b_is$numerator,
                                               den = tup_b_is$denominator,
                                               S = NULL,
                                               b = NULL,
                                               method = "mcmcse")
                     
                     out_b_tis <- elppd.details(num = tup_b_tis$numerator,
                                                den = tup_b_tis$denominator,
                                                S = NULL,
                                                b = NULL,
                                                method = "mcmcse")
                     
                     
                     out_b_psis <- elppd.details(num = tup_b_psis$numerator,
                                                 den = tup_b_psis$denominator,
                                                 S = NULL,
                                                 b = NULL,
                                                 method = "mcmcse")
                     
                     # Zanella outputs
                     
                     out_z_is <- elppd.details(num = tup_z_is$numerator,
                                               den = tup_z_is$denominator,
                                               S = NULL,
                                               b = NULL,
                                               method = "mcmcse")
                     
                     
                     out_z_tis <- elppd.details(num = tup_z_tis$numerator,
                                                den = tup_z_tis$denominator,
                                                S = NULL,
                                                b = NULL,
                                                method = "mcmcse")
                     
                     
                     out_z_psis <- elppd.details(num = tup_z_psis$numerator,
                                                 den = tup_z_psis$denominator,
                                                 S = NULL,
                                                 b = NULL,
                                                 method = "mcmcse")
                     
                     
                     c(elppd_g_is   = out_g_is$elppd,
                       var_g_is     = out_g_is$asymp_se,
                       
                       elppd_g_tis  = out_g_tis$elppd,
                       var_g_tis    = out_g_tis$asymp_se,
                       
                       elppd_g_psis = out_g_psis$elppd,
                       var_g_psis   = out_g_psis$asymp_se,
                       
                       acc_g = fit_g$acc,
                       
                       elppd_b_is   = out_b_is$elppd,
                       var_b_is     = out_b_is$asymp_se,
                       
                       elppd_b_tis  = out_b_tis$elppd,
                       var_b_tis    = out_b_tis$asymp_se,
                       
                       elppd_b_psis = out_b_psis$elppd,
                       var_b_psis   = out_b_psis$asymp_se,
                       
                       acc_b = fit_b$acc,
                       
                       elppd_z_is   = out_z_is$elppd,
                       var_z_is     = out_z_is$asymp_se,
                       
                       elppd_z_tis  = out_z_tis$elppd,
                       var_z_tis    = out_z_tis$asymp_se,
                       
                       elppd_z_psis = out_z_psis$elppd,
                       var_z_psis   = out_z_psis$asymp_se,
                       
                       acc_z = fit_z$acc
                     )
                     
                   }

stopCluster(cl)

results <- as.data.frame(results)

save(results,file = paste0(K,
                           "_fold_elppd_Gaussian.RData"))