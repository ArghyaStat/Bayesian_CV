# Stable log-sum-exp

log_sum_exp <- function(x){
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
}


# Truncated importance sampling

truncate_weights <- function(w){
  
  tau <- sqrt(length(w)) * mean(w)
  
  pmin(w, tau)
}

tuple_function <- function(log_lik_mat,
                           kappa,
                           method,
                           weight_type = "is",
                           alpha = NULL){
  
  method <- match.arg(method, c("Gelfand",
                                "Bronze",
                                "Zanella"))
  
  weight_type <- match.arg(weight_type, c("is",
                                          "tis",
                                          "psis"))
  
  S <- nrow(log_lik_mat)
  n <- ncol(log_lik_mat)
  
  if(length(kappa) != n){
    stop("length(kappa) must equal ncol(log_lik_mat)")
  }
  
  K <- max(kappa)

  
  fold_loglik <- matrix(0, S, K)
  
  for(k in seq_len(K)){
    
    idx <- which(kappa == k)
    
    fold_loglik[, k] <-
      rowSums(log_lik_mat[, idx, drop = FALSE])
  }
  
  
  total_loglik <- rowSums(log_lik_mat)
  
  
  # Raw importance weights
  
  w_mat <- matrix(0, S, n)
  
  # Gelfand
  
  
  if(method == "Gelfand"){
    
    for(i in seq_len(n)){
      
      k_i <- kappa[i]
      
      w_mat[, i] <-
        exp(-fold_loglik[, k_i])
    }
  }
  
  
  # Bronze
 
  if(method == "Bronze"){
    
    if(is.null(alpha)){
      stop("alpha must be supplied for Bronze.")
    }
    
    for(i in seq_len(n)){
      
      k_i <- kappa[i]
      
      log_test <-
        fold_loglik[, k_i]
      
      log_train <-
        total_loglik - log_test
      
      log_w <-
        (1 - alpha) * log_train -
        alpha * log_test
      
      w_mat[, i] <- exp(log_w)
    }
  }
  
 
  # Zanella
  
  if(method == "Zanella"){
    
    log_inv_fold <- -fold_loglik
    
    log_norm <-
      apply(log_inv_fold, 1, log_sum_exp)
    
    for(i in seq_len(n)){
      
      k_i <- kappa[i]
      
      w_mat[, i] <-
        exp(log_inv_fold[, k_i] -
              log_norm)
    }
  }
  
  # IS / TIS / PSIS
  
  for(i in seq_len(n)){
    
    if(weight_type == "is"){
      
      w <- w_mat[, i]
    }
    
    if(weight_type == "tis"){
      
      w <- truncate_weights(w_mat[, i])
    }
    
    if(weight_type == "psis"){
      
      psis_details <-
        loo::psis(log(w_mat[, i]),
                  r_eff = 1)
      
      w <- weights(psis_details,
                   log = FALSE,
                   normalize = FALSE)
    }
    
    w_mat[, i] <- w
  }
  
  # Tuple construction
  
  numerator <- exp(log_lik_mat) * w_mat
  denominator <- w_mat
  
  list(numerator = numerator,
       denominator = denominator)
}

# Batch means quadratic form t(v) %*% Omega %*% v

bm.quad.form <- function(Y, v,
                         b = NULL){
  
  Y <- as.matrix(Y)
  S <- nrow(Y)
  p <- ncol(Y)
  
  if(is.null(b)){ b <- floor(sqrt(S))}
  
  a <- floor(S / b)
  S_trim <- a * b
  
  Y <- Y[1:S_trim , , drop = FALSE]
  
  batch_means <- matrix(0, a, p)
  
  for(k in seq_len(a)){
    
    idx <- ((k - 1) * b + 1):(k * b)
    
    batch_means[k, ] <- colMeans( Y[idx, , drop = FALSE])
  }
  
  # Centered batch means matrix
  
  grand_mean <- colMeans(batch_means)
  
  B <- sqrt(b / (a - 1)) * sweep(batch_means, 2, grand_mean)
  
  Bv <- B %*% v
  
  as.numeric(crossprod(Bv))
}


# Asymptotic variance estimator

asymp.var.elppd <- function(num, den, S, b = NULL){
    
    v <- c(1 / colMeans(num), -1 / colMeans(den))
    Y <- cbind(num, den)
    mu <- colMeans(num) / colMeans(den)
    
  
  qf <- bm.quad.form(Y = Y, v = v, b = b)
  
  list(elppd = sum(log(mu)), asymp_var = qf / S)
}

# tuple_function <- function(log_lik_mat,
#                            method,
#                            weight_type = "is",
#                            alpha = NULL){
#   
#   method <- match.arg(method, c("Gelfand",
#                                 "Bronze",
#                                 "Zanella"))
#   
#   weight_type <- match.arg(weight_type, 
#                            c("is", "tis", "psis"))
#   
#   S <- nrow(log_lik_mat)
#   n <- ncol(log_lik_mat)
#   
#   # Construct raw weights
#   
#   w_mat <- matrix(0, S, n)
#   
#   # Gelfand
#   
#   if(method == "Gelfand"){
#     
#     w_mat <- exp(-log_lik_mat)
#   }
#   
#   # Bronze
#   
#   if(method == "Bronze"){
#     
#     if(is.null(alpha)){stop("alpha must be supplied for Bronze.")}
#     
#     w_mat <- exp(-alpha * log_lik_mat)
#   }
#   
#   # Zanella
#   
#   if(method == "Zanella"){
#     
#     log_inv <- -log_lik_mat
#     
#     log_norm <- apply(log_inv, 1, log_sum_exp)
#     
#     for(i in 1:n){
#       
#       w_mat[, i] <- exp(log_inv[, i] - log_norm)
#     }
#   }
#   
#   
#   for(i in 1:n){
#     
#     
#     # Ordinary IS
#     
#     if(weight_type == "is"){
#       
#       w <- w_mat[, i]
#     }
#     
#     
#     # Truncated IS
#     
#     if(weight_type == "tis"){
#       
#       w <- truncate_weights(w_mat[, i])
#     }
#     
#     # PSIS
#     
#     if(weight_type == "psis"){
#       
#       psis_details <- loo::psis(log(w_mat[, i]), r_eff = 1)
#       
#       w <- weights(psis_details, log=FALSE, normalize = FALSE)
#     }
#     
#     w_mat[, i] <- w
#   }
#   
#   # tuple construction
#   
#   numerator <- exp(log_lik_mat) * w_mat
#   denominator <- w_mat
#   
#   list(numerator = numerator,
#        denominator = denominator)
# }
