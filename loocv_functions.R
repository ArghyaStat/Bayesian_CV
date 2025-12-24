# Log posterior for simple linear regression model

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
  log_prior(theta) + sum(ll) + log(sum(exp(-ll))) # from supplement of Zanella's paper
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


psis_smooth <- function(log_ratios, tail_frac){
  S <- length(log_ratios)
  M <- ceiling(tail_frac * S)
  lw <- log_ratios - max(log_ratios)
  ord <- order(lw)
  tail <- ord[(S-M+1):S]
  y <- exp(lw[tail]) - min(exp(lw[tail]))
  k <- 0.5*(mean(y)^2/var(y) - 1)
  sigma <- mean(y)*(1+k)
  p <- (seq_len(M)-0.5)/M
  smoothed <- min(exp(lw[tail])) + (sigma/k)*((1-p)^(-k)-1)
  lw[tail] <- log(smoothed)
  lw <- pmin(lw, 0.5*log(S))
  w <- exp(lw)
  w/sum(w)
}
