# This file contains the importance sampling proposals of different methods
# "posterior", "Bronze", "mixture" 


# posterior proposal

log_q_post <- function(theta,
                    X,
                    y,
                    log_like,
                    log_prior) {
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + sum(ll)
}

grad_log_q_post <- function(theta,
                         X,
                         y,
                         grad_log_like,
                         grad_log_prior) {
  
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
                    log_prior) {
  
  ll <- log_like(theta, X, y)
  
  log_prior(theta) + alpha * sum(ll)
}

grad_log_q_b <- function(theta,
                         X,
                         y,
                         alpha,
                         grad_log_like,
                         grad_log_prior) {
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  grad_lp + alpha * colSums(grad_ll)
}

# mixture proposal

log_q_mix <- function(theta,
                    X,
                    y,
                    kappa,
                    log_like,
                    log_prior) {
  
  if (is.null(kappa) || length(kappa) == 0 || max(kappa) < 1) {
    stop("kappa must be a non-empty vector of positive integers (fold allocations).")
  }

  ll <- log_like(theta, X, y)
  
  K <- max(kappa)
  
  fold_loglik <- numeric(K)
  
  for(k in seq_len(K)) {
    
    idx <- which(kappa == k)
    
    fold_loglik[k] <- sum(ll[idx])
  }
  
  log_prior(theta) + sum(ll) + log_sum_exp(- fold_loglik)
}

grad_log_q_mix <- function(theta,
                         X,
                         y,
                         kappa,
                         grad_log_like,
                         grad_log_prior,
                         log_like) {
  
  if (is.null(kappa) || length(kappa) == 0 || max(kappa) < 1) {
    stop("kappa must be a non-empty vector of positive integers (fold allocations).")
  }

  ll <- log_like(theta, X, y)
  
  grad_ll <- grad_log_like(theta, X, y)
  
  grad_lp <- grad_log_prior(theta)
  
  K <- max(kappa)
  
  fold_loglik <- numeric(K)
  
  fold_grad <- matrix(0, nrow = K, ncol = ncol(grad_ll))
  
  for(k in seq_len(K)) {
    
    idx <- which(kappa == k)
    
    fold_loglik[k] <- sum(ll[idx])
    
    fold_grad[k, ] <- colSums(grad_ll[idx, ,drop = FALSE])
  }
  
  log_w <-  - fold_loglik
  
  log_w <- log_w - log_sum_exp(log_w)
  
  w <-  exp(log_w)
  
  grad_mix <- colSums(fold_grad * w)
  
  grad_lp + colSums(grad_ll) - grad_mix
}

## Potential energy functions

U_post <- function(theta,
                X,
                y,
                log_like,
                log_prior,
                ...) { - log_q_post(theta = theta,
           X = X,
           y = y,
           log_like = log_like,
           log_prior = log_prior)
}

grad_post <- function(theta,
                   X,
                   y,
                   grad_log_like,
                   grad_log_prior,
                   ...) { - grad_log_q_post(theta = theta,
                X = X,
                y = y,
                grad_log_like = grad_log_like,
                grad_log_prior = grad_log_prior)
}

U_b <- function(theta,
                X,
                y,
                alpha,
                log_like,
                log_prior,
                ...) { - log_q_b(theta = theta,
           X = X,
           y = y,
           alpha = alpha,
           log_like = log_like,
           log_prior = log_prior)
}

grad_b <- function(theta,
                   X,
                   y,
                   alpha,
                   grad_log_like,
                   grad_log_prior,
                   ...) { - grad_log_q_b(theta = theta,
                X = X,
                y = y,
                alpha = alpha,
                grad_log_like = grad_log_like,
                grad_log_prior = grad_log_prior)
}

U_mix <- function(theta,
                X,
                y,
                kappa,
                log_like,
                log_prior,
                ...) { - log_q_mix(theta = theta,
           X = X,
           y = y,
           kappa = kappa,
           log_like = log_like,
           log_prior = log_prior)
}

grad_mix <- function(theta,
                   X,
                   y,
                   kappa,
                   grad_log_like,
                   grad_log_prior,
                   log_like,
                   ...) { - grad_log_q_mix(theta = theta,
                X = X,
                y = y,
                kappa = kappa,
                grad_log_like = grad_log_like,
                grad_log_prior = grad_log_prior,
                log_like = log_like)
}