rm(list = ls())

library(MASS)
library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)



n <- 50
p <- 10
K <- 50
S <- 5e4
sigma2 <- 1
tau2 <- 1e2
reps <- 100


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




log_sum_exp <- function(x){
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
}


log_like <- function(theta, X, y){
  
  mu <- X %*% theta
  
  dnorm(y, mu, sqrt(sigma2), log = TRUE)
}



log_prior <- function(theta){
  
  sum(dnorm(theta, 0, sqrt(tau2), log = TRUE))
}


# Exact ELPPD

exact_elpd <- function(X, y){
  
  
  V_post <- solve(crossprod(X)/sigma2 + diag(1/tau2, p))
  
  m_post <- V_post %*% crossprod(X, y)/sigma2
  
  elpd <- 0
  
  for(i in 1:n){
    
    Xi <- matrix(X[i,], 1, p)
    
    X_minus <- X[-i,,drop = FALSE]
    y_minus <- y[-i]
    
    V_minus <- solve(
      crossprod(X_minus)/sigma2 +
        diag(1/tau2, p)
    )
    
    m_minus <- V_minus %*%
      crossprod(X_minus, y_minus)/sigma2
    
    pred_mean <- Xi %*% m_minus
    
    pred_var <- sigma2 +
      Xi %*% V_minus %*% t(Xi)
    
    elpd <- elpd +
      dnorm(
        y[i],
        pred_mean,
        sqrt(pred_var),
        log = TRUE
      )
  }
  
  as.numeric(elpd)
}

true_elpd <- exact_elpd(X, y)


# Gelfand

log_q_g <- function(theta){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll)
}

grad_log_q_g <- function(theta){
  
  mu <- X %*% theta
  
  r <- y - mu
  
  crossprod(X, r)/sigma2 -
    theta/tau2
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
  
  mu <- X %*% theta
  
  r <- y - mu
  
  (1 - 1/K) *
    crossprod(X, r)/sigma2 -
    theta/tau2
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
  
  mu <- X %*% theta
  
  r <- y - mu
  
  grad_ll <- crossprod(X, r)/sigma2
  
  ll_vec <- log_like(theta, X, y)
  
  log_w <- -ll_vec
  
  log_w <- log_w -
    log_sum_exp(log_w)
  
  w <- exp(log_w)
  
  grad_lse <- -crossprod(
    X,
    r * w
  )/sigma2
  
  grad_ll +
    grad_lse -
    theta/tau2
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




cl <- makeCluster(50)

registerDoParallel(cl)

results <- foreach(r = 1:reps, .combine = rbind, .packages = c("mcmcse", "MASS")) %dopar% {
  
  x_init <- lm(
    y ~ X - 1
  )$coef
  
  
  
  theta_g <- leapfrog_hmc(
    0.15, 20, S,
    U_g, grad_g,
    x_init, diag(p)
  )
  
  theta_b <- leapfrog_hmc(
    0.17, 20, S,
    U_b, grad_b,
    x_init, diag(p)
  )
  
  theta_z <- leapfrog_hmc(
    0.17, 20, S,
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


colnames(results) <- c(
  "elpd_g", "var_g", "acc_g",
  "elpd_b", "var_b", "acc_b",
  "elpd_z", "var_z", "acc_z"
)

save(
  results,
  true_elpd,
  file = paste0(
    K,
    "_fold_gaussian_exact_elpd.RData"
  )
)