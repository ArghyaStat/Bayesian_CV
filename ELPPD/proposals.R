# This file contains the importance sampling proposals of different methods
# "Gelfand", "Bronze", "Zanella" 


# Numerically stable log-sum-exp


log_sum_exp <- function(x){
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
}


# Gelfand proposal

log_q_g <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll)
}


grad_log_q_g <- function(theta,
                         X,
                         y,
                         grad_log_like,
                         grad_log_prior){
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  grad_lp + colSums(grad_ll)
}


# Bronze proposal

log_q_b <- function(theta,
                    X,
                    y,
                    alpha,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) +
    alpha * sum(ll)
}


grad_log_q_b <- function(theta,
                         X,
                         y,
                         alpha,
                         grad_log_like,
                         grad_log_prior){
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  grad_lp +
    alpha * colSums(grad_ll)
}



# Zanella proposal

log_q_z <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior){
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) +
    sum(ll) +
    log_sum_exp(-ll)
}


grad_log_q_z <- function(theta,
                         X,
                         y,
                         grad_log_like,
                         grad_log_prior,
                         log_like){
  
  ll <- log_like(theta, X, y)
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  log_w <- -ll
  
  log_w <- log_w - log_sum_exp(log_w)
  
  w <- exp(log_w)
  
  grad_mix <- colSums(
    grad_ll * w
  )
  
  grad_lp +
    colSums(grad_ll) -
    grad_mix
}



## Potential energy functions

U_g <- function(theta,
                X,
                y,
                log_like,
                log_prior){
  
  -log_q_g(
    theta = theta,
    X = X,
    y = y,
    log_like = log_like,
    log_prior = log_prior
  )
}


grad_g <- function(theta,
                   X,
                   y,
                   grad_log_like,
                   grad_log_prior){
  
  -grad_log_q_g(
    theta = theta,
    X = X,
    y = y,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior
  )
}



U_b <- function(theta,
                X,
                y,
                alpha,
                log_like,
                log_prior){
  
  -log_q_b(
    theta = theta,
    X = X,
    y = y,
    alpha = alpha,
    log_like = log_like,
    log_prior = log_prior
  )
}


grad_b <- function(theta,
                   X,
                   y,
                   alpha,
                   grad_log_like,
                   grad_log_prior){
  
  -grad_log_q_b(
    theta = theta,
    X = X,
    y = y,
    alpha = alpha,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior
  )
}



U_z <- function(theta,
                X,
                y,
                log_like,
                log_prior){
  
  -log_q_z(
    theta = theta,
    X = X,
    y = y,
    log_like = log_like,
    log_prior = log_prior
  )
}


grad_z <- function(theta,
                   X,
                   y,
                   grad_log_like,
                   grad_log_prior,
                   log_like){
  
  -grad_log_q_z(
    theta = theta,
    X = X,
    y = y,
    grad_log_like = grad_log_like,
    grad_log_prior = grad_log_prior,
    log_like = log_like
  )
}