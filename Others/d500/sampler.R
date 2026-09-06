## Leapfrog HMC for "Posterior" and "Bronze" IS

leapfrog_hmc <- function(epsilon,
                         L,
                         niters,
                         U,
                         grad_U,
                         x_init,
                         ...){
  
  d <- length(x_init)
  
  x <- matrix(NA, niters, d)
  
  x[1, ] <- x_init
  
  p <- matrix(rnorm(niters * d), niters, d)
  
  acc <- 0
  
  for(iter in 2:niters){
    
    x_prop <- x[iter - 1, ]
    
    p_prop <- p[iter, ]
    
    # Half momentum step
    
    g <- grad_U(x_prop, ...)
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    # Leapfrog steps
    
    for(l in 1:L){
      
      x_prop <- x_prop + epsilon * p_prop
      
      g <- grad_U(x_prop, ...)
      
      if(l != L){
        
        p_prop <- p_prop - epsilon * g
      }
    }
    
    # Final half momentum step
    
    p_prop <- p_prop - 0.5 * epsilon * g
    
    # Hamiltonians
    
    H_curr <- U(x[iter - 1, ], ...) + 0.5 * sum(p[iter, ]^2)
    
    H_prop <- U(x_prop, ... ) + 0.5 * sum(p_prop^2)
    
    # Accept/reject
    
    if(log(runif(1)) < -(H_prop - H_curr)){
      
      x[iter, ] <- x_prop
      acc <- acc + 1
      
    } else {
      
      x[iter, ] <- x[iter - 1, ]
    }
  }
  
  list(
    samples = x,
    acc = acc / niters
  )
}

non_rev_sampler <- function(niters,
                            epsilon,
                            L,
                            xi,
                            x_init,
                            M,
                            X = NULL,
                            y = NULL,
                            log_prior,
                            grad_log_prior,
                            log_like,
                            grad_log_like,
                            kappa,
                            modes = NULL,
                            R_chol = NULL,
                            log_det_V = NULL,
                            P_list = NULL,
                            B_list = NULL,
                            B_prime_list = NULL,
                            v_prior = NULL) {
  n <- if (is.null(y)) length(kappa) else length(y)
  d <- length(x_init)
  
  samples <- matrix(NA_real_, nrow = niters, ncol = d)
  samples[1, ] <- x_init
  
  I_samples <- integer(niters)
  I_samples[1] <- sample(1:n, 1)
  
  V <- matrix(sample(c(-1, 1), n * n, replace = TRUE), nrow = n, ncol = n)
  acc <- 0
  
  is_mat <- is.matrix(M)
  if (is_mat) {
    M_inv <- chol2inv(chol(M))
    chol_M <- chol(M)
  } else {
    M_inv <- 1 / M
    chol_M <- sqrt(M)
  }
  
  K <- max(kappa)
  fold_train_idx <- lapply(1:K, function(k) which(kappa != k))
  fold_test_idx <- lapply(1:K, function(k) which(kappa == k))
  
  # Pre-slice X and y for all K folds:
  X_tr_list <- lapply(fold_train_idx, function(idx) X[idx, , drop = FALSE])
  y_tr_list <- lapply(fold_train_idx, function(idx) y[idx])
  
  for (t in 2:niters) {
    theta_curr <- samples[t - 1, ]
    I_curr <- I_samples[t - 1]
    
    k_2 <- sample(n - 1, 1)
    if (k_2 >= I_curr) {
      k_2 <- k_2 + 1
    }
    
    idx_min <- min(I_curr, k_2)
    idx_max <- max(I_curr, k_2)
    v <- V[idx_min, idx_max]
    
    if (runif(1) < xi / n) {
      v <- -v
    }
    
    if (v == 1) {
      k_minus <- idx_min
      k_plus <- idx_max
    } else {
      k_minus <- idx_max
      k_plus <- idx_min
    }
    
    I_prop <- I_curr
    theta_prop <- theta_curr
    
    if (I_curr == k_plus) {
      v <- -v
    } else {
      fold_minus <- kappa[k_minus]
      fold_plus <- kappa[k_plus]
      
      if (fold_minus == fold_plus) {
        I_prop <- k_plus
      } else {
        if (!is.null(modes) && !is.null(P_list) && !is.null(B_list) && !is.null(B_prime_list) && !is.null(log_det_V) && !is.null(v_prior)) {
          z <- (1 / sqrt(v_prior)) * (theta_curr - modes[fold_minus, ])
          z <- z + P_list[[fold_minus]] %*% (B_prime_list[[fold_minus]] * as.vector(crossprod(P_list[[fold_minus]], z)))
          
          z_colored <- z + P_list[[fold_plus]] %*% (B_list[[fold_plus]] * as.vector(crossprod(P_list[[fold_plus]], z)))
          theta_prop <- modes[fold_plus, ] + sqrt(v_prior) * z_colored
          
          ll_curr <- log_like(theta_curr, X, y)
          ll_prop <- log_like(theta_prop, X, y)
          
          log_p_curr <- log_prior(theta_curr) + sum(ll_curr) - sum(ll_curr[fold_test_idx[[fold_minus]]])
          log_p_prop <- log_prior(theta_prop) + sum(ll_prop) - sum(ll_prop[fold_test_idx[[fold_plus]]])
          log_J <- 0.5 * log_det_V[fold_plus] - 0.5 * log_det_V[fold_minus]
          
          if (log(runif(1)) < log_p_prop + log_J - log_p_curr) {
            I_prop <- k_plus
            theta_curr <- theta_prop
          } else {
            v <- -v
            theta_prop <- theta_curr
          }
        } else if (!is.null(modes) && !is.null(R_chol) && !is.null(log_det_V)) {
          z <- solve(t(R_chol[[fold_minus]]), theta_curr - modes[fold_minus, ])
          theta_prop <- modes[fold_plus, ] + as.vector(t(R_chol[[fold_plus]]) %*% z)
          
          ll_curr <- log_like(theta_curr, X, y)
          ll_prop <- log_like(theta_prop, X, y)
          
          log_p_curr <- log_prior(theta_curr) + sum(ll_curr) - sum(ll_curr[fold_test_idx[[fold_minus]]])
          log_p_prop <- log_prior(theta_prop) + sum(ll_prop) - sum(ll_prop[fold_test_idx[[fold_plus]]])
          log_J <- 0.5 * log_det_V[fold_plus] - 0.5 * log_det_V[fold_minus]
          
          if (log(runif(1)) < log_p_prop + log_J - log_p_curr) {
            I_prop <- k_plus
            theta_curr <- theta_prop
          } else {
            v <- -v
            theta_prop <- theta_curr
          }
        } else {
          ll_all <- log_like(theta_curr, X, y)
          log_p_plus <- sum(ll_all) - sum(ll_all[fold_test_idx[[fold_plus]]])
          log_p_minus <- sum(ll_all) - sum(ll_all[fold_test_idx[[fold_minus]]])
          
          if (log(runif(1)) < log_p_plus - log_p_minus) {
            I_prop <- k_plus
          } else {
            v <- -v
          }
        }
      }
    }
    
    if (runif(1) < xi / n) {
      v <- -v
    }
    V[idx_min, idx_max] <- v
    I_curr <- I_prop
    I_samples[t] <- I_curr
    
    active_fold <- kappa[I_curr]
    X_tr <- X_tr_list[[active_fold]]
    y_tr <- y_tr_list[[active_fold]]
    
    if (is_mat) {
      p_curr <- as.vector(rnorm(d) %*% chol_M)
    } else {
      p_curr <- rnorm(d) * chol_M
    }
    
    theta_prop <- theta_curr
    p_prop <- p_curr
    
    g <- grad_U_fold(theta_prop, X_tr, y_tr)
    p_prop <- p_prop - 0.5 * epsilon * g
    
    for (l in 1:L) {
      if (is_mat) {
        theta_prop <- theta_prop + epsilon * as.vector(M_inv %*% p_prop)
      } else {
        theta_prop <- theta_prop + epsilon * (M_inv * p_prop)
      }
      
      g <- grad_U_fold(theta_prop, X_tr, y_tr)
      if (l != L) {
        p_prop <- p_prop - epsilon * g
      }
    }
    p_prop <- p_prop - 0.5 * epsilon * g
    
    if (is_mat) {
      kin_curr <- 0.5 * sum(p_curr * as.vector(M_inv %*% p_curr))
      kin_prop <- 0.5 * sum(p_prop * as.vector(M_inv %*% p_prop))
    } else {
      kin_curr <- 0.5 * sum(p_curr^2 * M_inv)
      kin_prop <- 0.5 * sum(p_prop^2 * M_inv)
    }
    
    H_curr <- U_fold(theta_curr, X_tr, y_tr) + kin_curr
    H_prop <- U_fold(theta_prop, X_tr, y_tr) + kin_prop
    
    if (log(runif(1)) < -(H_prop - H_curr)) {
      samples[t, ] <- theta_prop
      acc <- acc + 1
    } else {
      samples[t, ] <- theta_curr
    }
  }
  
  list(
    samples = samples,
    I_samples = I_samples,
    acc = acc / (niters - 1)
  )
}
