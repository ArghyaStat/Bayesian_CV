rm(list = ls())

library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)

# Simulation design

n <- 50
p <- 10
K <- 10
S <- 1e4
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

# K-fold splitting

fold_id <- sample(rep(1:K, length.out = n))

# indexing map kappa(i)

kappa <- function(i, fold_id){
  fold_id[i]
}

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

# Gradients for different methods

grad_post <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  -( t(X) %*% (y - p_hat) + theta/sigma2 )
}

grad_b <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  -( (K - 1)/K * t(X) %*% (y - p_hat) + theta/sigma2 )
}

grad_z <- function(theta){
  
  eta <- X %*% theta
  
  p_hat <- 1/(1 + exp(-eta))
  
  r <- y - p_hat
  
  grad_ll <- t(X) %*% r
  
  ll_vec <- log_like(theta, X, y)
  
  fold_ll <- sapply(1:K, function(k){
    
    idx <- which(fold_id == k)
    
    sum(ll_vec[idx])
  })
  
  log_w <- -fold_ll
  
  log_w <- log_w - log_sum_exp(log_w)
  
  w <- exp(log_w)
  
  grad_mix <- 0
  
  for(k in 1:K){
    
    idx <- which(fold_id == k)
    
    grad_k <- t(X[idx,,drop = FALSE]) %*% r[idx]
    
    grad_mix <- grad_mix - w[k] * grad_k
  }
  
  -(grad_ll + grad_mix + theta/sigma2)
}

# Potentials

U_post <- function(theta){
  
  -(sum(log_like(theta, X, y)) + log_prior(theta))
}

U_b <- function(theta){
  
  ll <- sum(log_like(theta, X, y))
  
  -( (K - 1)/K * ll + log_prior(theta) )
}

U_z <- function(theta){
  
  ll_vec <- log_like(theta, X, y)
  
  fold_ll <- sapply(1:K, function(k){
    
    idx <- which(fold_id == k)
    
    sum(ll_vec[idx])
  })
  
  -( sum(ll_vec) + log_sum_exp(-fold_ll) + log_prior(theta) )
}

# HMC code from (HMC for dummies)

leapfrog_hmc <- function(epsilon,
                         L,
                         niters,
                         U,
                         grad,
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
    
    g <- grad(x_prop)
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    for(l in 1:L){
      
      x_prop <- x_prop + epsilon * (M_inv %*% p_prop)
      
      g <- grad(x_prop)
      
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

# fold likelihoods

compute_fold_lik <- function(lik){
  
  sapply(1:K, function(k){
    
    idx <- which(fold_id == k)
    
    apply(lik[,idx,drop = FALSE], 1, prod)
  })
}


reps <- 100

cl <- makeCluster(50)

registerDoParallel(cl)


results <- foreach(r = 1:reps,
                   .combine = rbind,
                   .packages = c("mcmcse")) %dopar% {
                     
                     x_init <- glm(y ~ X - 1,
                                   family = "binomial")$coef
                     
                     # -------------------------
                     # HMC
                     # -------------------------
                     theta_g <- leapfrog_hmc(
                       0.35, 10, S,
                       U_post, grad_post,
                       x_init, diag(p)
                     )
                     
                     theta_b <- leapfrog_hmc(
                       0.35, 10, S,
                       U_b, grad_b,
                       x_init, diag(p)
                     )
                     
                     theta_z <- leapfrog_hmc(
                       0.15, 10, S,
                       U_z, grad_z,
                       x_init, diag(p)
                     )
                     
                     acc_g <- theta_g$acc
                     acc_b <- theta_b$acc
                     acc_z <- theta_z$acc
                     
                     samples_g <- theta_g$samples
                     samples_b <- theta_b$samples
                     samples_z <- theta_z$samples
                     
                     # Likelihood computations from different methods
                     
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
                     
                     fold_lik_b <- compute_fold_lik(lik_b)
                     
                     fold_lik_z <- compute_fold_lik(lik_z)
                     
                     full_lik_b <- apply(lik_b, 1, prod)
                     
                     # Gelfand's method
                     
                     mu_g <- numeric(n)
                     
                     num_g <- matrix(1, S, n)
                     
                     den_g <- 1/lik_g
                     
                     for(i in 1:n){
                       
                       mu_g[i] <- 1/mean(den_g[,i])
                     }
                     
                     elpd_g <- sum(log(mu_g))
                     
                     grad_g <- c(rep(1,n), -mu_g)
                     
                     Omega_g <- mcse.multi(
                       cbind(num_g, den_g),
                       r = 1,
                       adjust = FALSE,
                       size = 100
                     )$cov
                     
                     var_g <- as.numeric(
                       t(grad_g) %*% Omega_g %*% grad_g / S
                     )
                     
                     # Bronze method
                     
                     mu_b <- numeric(n)
                     
                     num_b <- matrix(NA, S, n)
                     
                     den_b <- matrix(NA, S, n)
                     
                     for(i in 1:n){
                       
                       k_i <- kappa(i, fold_id)
                       
                       lik_fold <- fold_lik_b[,k_i]
                       
                       lik_train <- full_lik_b / lik_fold
                       
                       w <- (lik_train)^(1/K) *
                         (lik_fold)^(-(K - 1)/K)
                       
                       num_b[,i] <- lik_b[,i] * w
                       
                       den_b[,i] <- w
                       
                       mu_b[i] <- sum(num_b[,i]) /
                         sum(den_b[,i])
                     }
                     
                     elpd_b <- sum(log(mu_b))
                     
                     grad_b_vec <- c(
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
                       t(grad_b_vec) %*% Omega_b %*% grad_b_vec / S
                     )
                     
                    # Zanella's method
                     
                     mu_z <- numeric(n)
                     
                     num_z <- matrix(NA, S, n)
                     
                     den_z <- matrix(NA, S, n)
                     
                     denom_all <- rowSums(1/fold_lik_z)
                     
                     for(i in 1:n){
                       
                       k_i <- kappa(i, fold_id)
                       
                       w <- (1/fold_lik_z[,k_i]) /
                         denom_all
                       
                       num_z[,i] <- lik_z[,i] * w
                       
                       den_z[,i] <- w
                       
                       mu_z[i] <- sum(num_z[,i]) /
                         sum(den_z[,i])
                     }
                     
                     elpd_z <- sum(log(mu_z))
                     
                     grad_z_vec <- c(
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
                       t(grad_z_vec) %*% Omega_z %*% grad_z_vec / S
                     )
                     
                     c(elpd_g, var_g, acc_g,
                       elpd_b, var_b, acc_b,
                       elpd_z, var_z, acc_z)
                   }

stopCluster(cl)

# Outputs

colnames(results) <- c(
  "elpd_g", "var_g", "acc_g",
  "elpd_b", "var_b", "acc_b",
  "elpd_z", "var_z", "acc_z"
)

save(
  results,
  file = paste0(K, "_fold_elppd_with_var.RData")
)