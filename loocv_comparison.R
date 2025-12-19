rm(list = ls())
set.seed(1)

library(mcmcse)

#  Data Simulation

n <- 50
x <- rnorm(n)
y <- 1 + 2*x + rnorm(n)
sigma2 <- 1
tau2 <- 100

# Log posterior

log_prior <- function(theta){
  sum(dnorm(theta, 0, sqrt(tau2), log = TRUE))
}

log_like <- function(theta){
  mu <- theta[1] + theta[2]*x
  dnorm(y, mu, sqrt(sigma2), log = TRUE)
}

log_posterior <- function(theta){
  log_prior(theta) + sum(log_like(theta))
}

# potential and gradient (Gelfand)

log_post <- function(theta){
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- beta0 + beta1 * x
  
  ll <- sum(dnorm(y, mu, sqrt(sigma2), log = TRUE))
  lp <- sum(dnorm(theta, 0, sqrt(tau2), log = TRUE))
  
  ll + lp
}

U_post <- function(theta){
  -log_post(theta)
}

grad_U_post <- function(theta){
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu <- beta0 + beta1 * x
  r <- y - mu
  
  # gradient of log-likelihood
  grad_ll <- c(sum(r), sum(r * x)) / sigma2
  
  # gradient of log-prior
  grad_lp <- -theta / tau2
  
  -(grad_ll + grad_lp)
}



# Zanella's q_mix and gradient

log_q_mix <- function(theta){
 
  ll <- log_like(theta)
  log_prior(theta) + sum(ll) + log(sum(exp(-ll))) # from supplemental of paper
}

grad_log_q_mix <- function(theta){
  
  mu <- theta[1] + theta[2]*x
  r <- y - mu # residual
  
  # gradient of sum log-likelihood
  grad_ll <- c(sum(r), sum(r*x)) / sigma2
  
  ll <- log_like(theta)
  w <- exp(-ll) / sum(exp(-ll))
  
  grad_lse <- -colSums(cbind(r, r*x) * w) / sigma2
  grad_ll + grad_lse - theta/tau2
}

U_mix <- function(theta) -log_q_mix(theta)
grad_U_mix <- function(theta) -grad_log_q_mix(theta)

# HMC sampler (copied from our paper on HMC for dummies)

leapfrog_hmc <- function(epsilon, L, niters, U, grad, x_init, M) {
  d <- length(x_init)
  M_inv <- chol2inv(chol(M))
  x <- matrix(NA, niters, d)
  x[1,] <- x_init
  
  chol_M <- chol(M)
  p <- matrix(rnorm(niters*d), niters, d) %*% chol_M
  acc <- 0
  
  for(k in 2:niters){
    x_prop <- x[k-1,]
    p_prop <- p[k,]
    
    g <- grad(x_prop)
    p_prop <- p_prop - 0.5*epsilon*g
    
    for(l in 1:L){
      x_prop <- x_prop + epsilon*(M_inv %*% p_prop)
      g <- grad(x_prop)
      if(l != L) p_prop <- p_prop - epsilon*g
    }
    p_prop <- p_prop - 0.5*epsilon*g
    
    H_curr <- U(x[k-1,]) + 0.5*t(p[k,]) %*% (M_inv %*% p[k,])
    H_prop <- U(x_prop) + 0.5*t(p_prop) %*% (M_inv %*% p_prop)
    
    if(log(runif(1)) < -(H_prop - H_curr)){
      x[k,] <- x_prop
      acc <- acc + 1
    } else {
      x[k,] <- x[k-1,]
    }
  }
  list(samples = x, acc = acc/niters)
}


# Gelfand sampler via HMC 

S <- 5e3

hmc_gelfand <- leapfrog_hmc(
  epsilon = 0.25,
  L = 20,
  niters = S,
  U = U_post,
  grad = grad_U_post,
  x_init = c(1, 2),
  M = diag(2)
)

acc_gelfand <- hmc_gelfand$acc
acc_gelfand

theta_g <- hmc_gelfand$samples


# Zanella (HMC on q_mix)
hmc_zanella <- leapfrog_hmc(
  epsilon = 0.24, L = 20, niters = S,
  U = U_mix, grad = grad_U_mix,
  x_init = c(1,2), M = diag(2)
)

acc_zanella <- hmc_zanella$acc
acc_zanella
theta_z <- hmc_zanella$samples


# Log-likelihood matrices (FIXED) of dimension S \times n


loglik_g <- loglik_z <- matrix(NA, S, n)

for(s in 1:S){
  loglik_g[s,] <- log_like(theta_g[s,])
  loglik_z[s,] <- log_like(theta_z[s,])
}

lik_g <- exp(loglik_g)
lik_z <- exp(loglik_z)

# LOOCV estimators

# mu_hat and ELPD from Gelfand ----

mu_g <- sapply(1:n, function(i) 1/mean(1/lik_g[,i]))
elpd_g <- sum(log(mu_g))


## Important sampling SE

# Gelfand
num_g <- matrix(1, S, n)
den_g <- 1/lik_g
grad_g <- c(rep(1, n), -mu_g)

vec_g <- cbind(num_g, den_g)

foo_g <- mcse.multi(vec_g, r = 1, adjust = FALSE, size = 100)
Omega_hat_g <- foo_g$cov
se_g <- as.numeric(t(grad_g) %*% Omega_hat_g %*% grad_g / S)

# mu_hat and ELPD from Zanella (self-normalized mixture) ----

mu_z <- numeric(n)
num_z <- den_z <- matrix(NA, S, n)

for(i in 1:n){
  w <- exp(-loglik_z[,i])
  w <- w / sum(w)
  num_z[,i] <- lik_z[,i] * w
  den_z[,i] <- w
  mu_z[i] <- sum(num_z[,i])/sum(w)
}

elpd_z <- sum(log(mu_z))


# Zanella (self-normalized)
g1 <- 1/colMeans(num_z)
g2 <- -1/colMeans(den_z)
grad_z <- c(g1, g2)

foo_z <- mcse.multi(cbind(num_z, den_z), r = 1, adjust = FALSE, size = 100)
Omega_hat_z <- foo_z$cov
se_z <- as.numeric(t(grad_z) %*% Omega_hat_z %*% grad_z / S)


cat("\nGelfand ELPD:", elpd_g, " SE:", se_g, "\n")
cat("Zanella ELPD:", elpd_z, " SE:", se_z, "\n")


reps <- 1e2

elpd_rep_g <- numeric(reps)
elpd_rep_z <- numeric(reps)

for (r in 1:reps) {
  
  # cat("Replication:", r, "\n")
  
  # Gelfand
  
  theta_g <- leapfrog_hmc(
    epsilon = 0.25,
    L = 20,
    niters = S,
    U = U_post,
    grad = grad_U_post,
    x_init = c(1, 2),
    M = diag(2)
  )$samples
  
  lik_g <- matrix(NA, S, n)
  for (s in 1:S) {
    lik_g[s, ] <- exp(log_like(theta_g[s, ]))
  }
  
  loocv_g <- numeric(n)
  for (i in 1:n) {
    loocv_g[i] <- 1 / mean(1 / lik_g[, i])
  }
  
  elpd_rep_g[r] <- sum(log(loocv_g))
  
  
  # Zanella (q_mix via HMC)
 
  
  theta_z <- leapfrog_hmc(
    epsilon = 0.24, L = 20, niters = S,
    U = U_mix, grad = grad_U_mix,
    x_init = c(1,2), M = diag(2)
  )$samples
  
  lik_z <- matrix(NA, S, n)
  for (s in 1:S) {
    lik_z[s, ] <- exp(log_like(theta_z[s, ]))
  }
  
  loocv_z <- numeric(n)
  for (i in 1:n) {
    w <- 1 / lik_z[, i]
    w <- w / sum(w)   # self-normalized mixture weights
    loocv_z[i] <- sum(lik_z[, i] * w)/sum(w)
  }
  
  elpd_rep_z[r] <- sum(log(loocv_z))
}



cat("\nReplication results:\n\n")

cat("Gelfand:\n")
cat("Mean ELPD:", mean(elpd_rep_g), "\n")
cat("Var  ELPD:", var(elpd_rep_g), "\n\n")

cat("Zanella:\n")
cat("Mean ELPD:", mean(elpd_rep_z), "\n")
cat("Var  ELPD:", var(elpd_rep_z), "\n")