# Mixture proposal log-density

log_q_mix <- function(theta, X, y, kappa, log_like, log_prior) {
  
  if (is.null(kappa) || length(kappa) == 0 || max(kappa) < 1) {
    stop("kappa must be a non-empty vector of positive integers (fold allocations).")
  }
  
  ll  <- log_like(theta, X, y)
  fold_loglik <- as.vector(rowsum.default(ll, kappa, reorder = FALSE))
  log_prior(theta) + sum(ll) + log_sum_exp(- fold_loglik)
}

# Mixture proposal gradient
grad_log_q_mix <- function(theta, X, y, kappa, grad_log_like, grad_log_prior, log_like) {
  
  if (is.null(kappa) || length(kappa) == 0 || max(kappa) < 1) {
    stop("kappa must be a non-empty vector of positive integers (fold allocations).")
  }
  
  ll  <- log_like(theta, X, y)
  grad_ll <- grad_log_like(theta, X, y)
  grad_lp <- grad_log_prior(theta)
  
  fold_loglik <- as.vector(rowsum.default(ll, kappa, reorder = FALSE))
  fold_grad <- rowsum.default(grad_ll, kappa, reorder = FALSE)
  
  log_w <- - fold_loglik
  log_w <- log_w - log_sum_exp(log_w)
  w <- exp(log_w)
  
  grad_mix <- colSums(fold_grad * w)
  grad_lp + colSums(grad_ll) - grad_mix
}

# Potential energy function and gradient for HMC leapfrog
U_mix <- function(theta, X, y, kappa, log_like, log_prior, ...) {
  
  - log_q_mix(theta = theta, X = X, y = y, kappa = kappa,
              log_like = log_like, log_prior = log_prior)
}

grad_mix <- function(theta, X, y, kappa, grad_log_like, grad_log_prior, log_like, ...) {
  
  - grad_log_q_mix(theta = theta, X = X, y = y, kappa = kappa,
                   grad_log_like = grad_log_like, 
                   grad_log_prior = grad_log_prior,
                   log_like = log_like)
}


# Leapfrog Hamiltonian Monte Carlo sampler
leapfrog_hmc <- function(epsilon, L, niters, U, grad_U, x_init, ...) {
  
  d    <- length(x_init)
  x    <- matrix(NA_real_, niters, d)
  x[1, ] <- x_init
  p    <- matrix(rnorm(niters * d), niters, d)
  acc  <- 0
  
  for (iter in 2:niters) {
    x_prop <- x[iter - 1, ]
    p_prop <- p[iter, ]
    
    g      <- grad_U(x_prop, ...)
    p_prop <- p_prop - 0.5 * epsilon * g
    
    for (l in 1:L) {
      x_prop <- x_prop + epsilon * p_prop
      g      <- grad_U(x_prop, ...)
      if (l != L) {
        p_prop <- p_prop - epsilon * g
      }
    }
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    H_curr <- U(x[iter - 1, ], ...) + 0.5 * sum(p[iter, ]^2)
    H_prop <- U(x_prop, ...)       + 0.5 * sum(p_prop^2)
    
    if (log(runif(1)) < -(H_prop - H_curr)) {
      x[iter, ] <- x_prop
      acc        <- acc + 1
    } else {
      x[iter, ] <- x[iter - 1, ]
    }
  }
  
  list(samples = x, acc = acc / (niters - 1))
}

# Observation-level log-likelihood and gradient
log_like <- function(theta, X, y) {
  mu <- X %*% theta
  dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE)
}

grad_log_like <- function(theta, X, y) {
  mu <- X %*% theta
  r  <- y - mu
  X * as.vector(r) / sigma2
}

log_prior <- function(theta) {
  sum(dnorm(theta, mean = 0, sd = sqrt(tau2), log = TRUE))
}

grad_log_prior <- function(theta) {
  -theta / tau2
}


#####
# Run now
d      <- 250
n      <- 500
eta    <- 0.1
sigma2 <- 1
tau2   <- 100
K     <- n
kappa <- 1:n

# HMC parameters
eps_mix <- 0.48
L_steps <- 20

X_curr      <- sim_dat$X
dropped_col <- sim_dat$dropped_col
y           <- sim_dat$y

fit_lm    <- lm(y ~ X_curr - 1)
beta_init <- as.vector(coef(fit_lm))
beta_init[is.na(beta_init)] <- 0
S_max  <- 5e3

fit_mix    <- leapfrog_hmc(
  epsilon = eps_mix, L = L_steps, niters = S_max,
  U = U_mix, grad_U = grad_mix, x_init = beta_init,
  X = X_curr, y = y, kappa = kappa,
  log_prior = log_prior,
  grad_log_prior = grad_log_prior,
  log_like = log_like,
  grad_log_like = grad_log_like)


library(Rcpp)
library(RcppArmadillo)
sourceCpp("sampler.cpp")
# Run now
d      <- 250
n      <- 500
eta    <- 0.1
sigma2 <- 1
tau2   <- 100
K      <- n
kappa  <- 1:n

# HMC parameters
eps_mix <- 0.48
L_steps <- 20

X_curr      <- sim_dat$X
dropped_col <- sim_dat$dropped_col
y           <- sim_dat$y

fit_lm    <- lm(y ~ X_curr - 1)
beta_init <- as.vector(coef(fit_lm))
beta_init[is.na(beta_init)] <- 0

S_max <- 1e4

system.time(fit_mix2 <- hmc_loocv_cpp(
  X = X_curr,
  y = y,
  theta_init = beta_init,
  sigma2 = sigma2,
  tau2 = tau2,
  epsilon = eps_mix,
  L = 10,
  niters = S_max
))


plot(acf(fit_mix$samples[,1]), type = "n")
for(i in 1:dim(fit_mix$samples)[2])
{
  foo <- acf(fit_mix$samples[,i], plot = FALSE)$acf
  lines(foo)
}

i <- 1
i <- i+1
plot(fit_mix2$samples[,i], fit_mix2$samples[,i+1], 
     pch = 16, col = "blue")
for(i in 11:31)
{
  plot(density(fit_mix2$samples[,i]))
}



