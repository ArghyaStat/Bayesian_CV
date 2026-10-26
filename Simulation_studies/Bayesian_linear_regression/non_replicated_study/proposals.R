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
