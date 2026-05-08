rm(list = ls())

library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)

# Simulation design

n <- 50
p <- 10
K <- 5
S <- 5e4
sigma2 <- 1e2

# Data generation

sim_data <- function(n, p){
  
  X <- matrix(rnorm(n*p), n, p)
  
  beta <- rnorm(p)
  
  eta <- X %*% beta
  
  prob <- 1/(1 + exp(-eta))
  
  y <- rbinom(n, 1, prob)
  
  list(X = X, y = y)
}

dat <- sim_data(n, p)

X <- dat$X
y <- dat$y

# Stable computation at log-scale adjustment

log_sum_exp <- function(x){
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
}

log_like <- function(theta, X, y){
  
  eta <- X %*% theta
  
  y * eta - log(1 + exp(eta))
}

log_prior <- function(theta){
  
  sum(dnorm(theta, 0, sqrt(sigma2), log = TRUE))
}


# Gelfand


log_q_g <- function(theta){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll)
}

grad_log_q_g <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  r <- y - p_hat
  
  t(X) %*% r - theta/sigma2
}

U_g <- function(theta){
  
  -log_q_g(theta)
}

grad_g <- function(theta){
  
  -grad_log_q_g(theta)
}


# Bronze


log_q_b <- function(theta){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) +
    (1 - 1/K) * sum(ll)
}

grad_log_q_b <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  r <- y - p_hat
  
  (1 - 1/K) * t(X) %*% r -
    theta/sigma2
}

U_b <- function(theta){
  
  -log_q_b(theta)
}

grad_b <- function(theta){
  
  -grad_log_q_b(theta)
}


# Zanella

log_q_mix <- function(theta){
  
  ll_vec <- log_like(theta, X, y)
  
  log_prior(theta) +
    sum(ll_vec) +
    log_sum_exp(-ll_vec)
}

grad_log_q_mix <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  r <- y - p_hat
  
 
  grad_ll <- t(X) %*% r
  
  ll_vec <- log_like(theta, X, y)
  
  # normalized weights
  log_w <- -ll_vec
  
  log_w <- log_w - log_sum_exp(log_w)
  
  w <- exp(log_w)
  
  # gradient of log-sum-exp correction
  grad_lse <- -colSums(
    X * as.vector(r * w)
  )
  
  grad_ll + grad_lse - theta/sigma2
}

U_z <- function(theta){
  
  -log_q_mix(theta)
}

grad_z <- function(theta){
  
  -grad_log_q_mix(theta)
}


# HMC

leapfrog_hmc <- function(epsilon,
                         L,
                         niters,
                         U,
                         grad_U,
                         x_init,
                         M){
  
  d <- length(x_init)
  
  M_inv <- chol2inv(chol(M))
  
  x <- matrix(NA, niters, d)
  
  x[1,] <- x_init
  
  chol_M <- chol(M)
  
  p <- matrix(rnorm(niters * d), niters, d) %*% chol_M
  
  acc <- 0
  
  for(iter in 2:niters){
    
    x_prop <- x[iter - 1,]
    
    p_prop <- p[iter,]
    
    g <- grad_U(x_prop)
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    for(l in 1:L){
      
      x_prop <- x_prop +
        epsilon * (M_inv %*% p_prop)
      
      g <- grad_U(x_prop)
      
      if(l != L){
        p_prop <- p_prop - epsilon * g
      }
    }
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    H_curr <- U(x[iter - 1,]) +
      0.5 * t(p[iter,]) %*% (M_inv %*% p[iter,])
    
    H_prop <- U(x_prop) +
      0.5 * t(p_prop) %*% (M_inv %*% p_prop)
    
    if(log(runif(1)) < -(H_prop - H_curr)){
      
      x[iter,] <- x_prop
      
      acc <- acc + 1
      
    } else {
      
      x[iter,] <- x[iter - 1,]
    }
  }
  
  list(samples = x,
       acc = acc/niters)
}


reps <- 100

cl <- makeCluster(50)

registerDoParallel(cl)

results <- foreach(r = 1:reps, .combine = rbind, .packages = c("mcmcse")) %dopar% {
  
          x_init <- glm(
            y ~ X - 1,
            family = "binomial"
          )$coef
          
         
          
          theta_g <- leapfrog_hmc(
            0.6, 10, S,
            U_g, grad_g,
            x_init, diag(p)
          )
          
          theta_b <- leapfrog_hmc(
            0.6, 10, S,
            U_b, grad_b,
            x_init, diag(p)
          )
          
          theta_z <- leapfrog_hmc(
            0.85, 10, S,
            U_z, grad_z,
            x_init, diag(p)
          )
          
          acc_g <- theta_g$acc
          acc_b <- theta_b$acc
          acc_z <- theta_z$acc
          
          samples_g <- theta_g$samples
          samples_b <- theta_b$samples
          samples_z <- theta_z$samples
          
          
          # Likelihoods matrices of dimension S \times n
          
          lik_g <- t(
            exp(
              sapply(1:S, function(s)
                log_like(samples_g[s,], X, y))
            )
          )
          
          lik_b <- t(
            exp(
              sapply(1:S, function(s)
                log_like(samples_b[s,], X, y))
            )
          )
          
          lik_z <- t(
            exp(
              sapply(1:S, function(s)
                log_like(samples_z[s,], X, y))
            )
          )
          
          
          # Gelfand
         
          
          mu_g <- numeric(n)
          
          den_g <- 1/lik_g
          
          for(i in 1:n){
            
            mu_g[i] <- 1/mean(den_g[,i])
          }
          
          elpd_g <- sum(log(mu_g))
          
          v_g <- -1/colMeans(den_g)
          
          Omega_g <- mcse.multi(
            den_g,
            r = 1,
            adjust = FALSE,
            size = 100
          )$cov
          
          var_g <- as.numeric(
            t(v_g) %*% Omega_g %*% v_g / S
          )
          
          
          # Bronze
          
          
          mu_b <- numeric(n)
          
          num_b <- lik_b^(1/K)
          
          den_b <- lik_b^(-(K - 1)/K)
          
          for(i in 1:n){
            
            mu_b[i] <- sum(num_b[,i]) /
              sum(den_b[,i])
          }
          
          elpd_b <- sum(log(mu_b))
          
          v_b <- c(
            1/colMeans(num_b),
            -1/colMeans(den_b)
          )
          
          Omega_b <- mcse.multi(
            cbind(num_b, den_b),
            r = 1,
            adjust = FALSE,
            size = 100
          )$cov
          
          var_b <- as.numeric(
            t(v_b) %*% Omega_b %*% v_b / S
          )
          
         
          # Zanella
          
          
          mu_z <- numeric(n)
          
          inv_lik <- 1/lik_z
          
          denom_all <- rowSums(inv_lik)
          
          num_z <- lik_z * (inv_lik / denom_all)
          
          den_z <- inv_lik / denom_all
          
          for(i in 1:n){
            
            mu_z[i] <- sum(num_z[,i]) /
              sum(den_z[,i])
          }
          
          elpd_z <- sum(log(mu_z))
          
          v_z <- c(
            1/colMeans(num_z),
            -1/colMeans(den_z)
          )
          
          Omega_z <- mcse.multi(
            cbind(num_z, den_z),
            r = 1,
            adjust = FALSE,
            size = 100
          )$cov
          
          var_z <- as.numeric(
            t(v_z) %*% Omega_z %*% v_z / S
          )
          
          c(
            elpd_g, var_g, acc_g,
            elpd_b, var_b, acc_b,
            elpd_z, var_z, acc_z
          )
}

stopCluster(cl)

# =========================================================
# OUTPUTS
# =========================================================

colnames(results) <- c(
  "elpd_g", "var_g", "acc_g",
  "elpd_b", "var_b", "acc_b",
  "elpd_z", "var_z", "acc_z"
)

save(
  results,
  file = paste0(K, "_fold_elppd_outputs.RData")
)