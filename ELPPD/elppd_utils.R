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

gpdfitnew <- function(x,
                      sort = TRUE){
  
  if(!is.vector(x) || length(x) <= 1)
    stop("Invalid input array.")
  
  ##########################################
  # Sorting
  ##########################################
  
  if(isTRUE(sort)){
    
    ord <- order(x)
    xsorted <- FALSE
    
  } else if(isFALSE(sort)) {
    
    xsorted <- TRUE
    
  } else {
    
    ord <- sort
    xsorted <- FALSE
  }
  
  n <- length(x)
  
  PRIOR <- 3
  
  m <- 30 + floor(sqrt(n))
  
  ##########################################
  # Candidate b values
  ##########################################
  
  bs <- 1:m
  
  bs <- 1 - sqrt(m / (bs - 0.5))
  
  if(xsorted){
    
    bs <-
      bs /
      (PRIOR *
       x[floor(n/4 + 0.5)]) +
      1 / x[n]
    
  } else {
    
    bs <-
      bs /
      (PRIOR *
       x[ord[floor(n/4 + 0.5)]]) +
      1 / x[ord[n]]
  }
  
  ##########################################
  # Compute k estimates
  ##########################################
  
  ks <- numeric(length(bs))
  
  for(j in seq_along(bs)){
    
    temp <- log1p(
      -bs[j] * x
    )
    
    ks[j] <- mean(temp)
  }
  
  ##########################################
  # Marginal posterior weights
  ##########################################
  
  L <- -n * (
    log(-bs / ks) +
    ks + 1
  )
  
  temp <- outer(L, L, "-")
  
  w <- 1 / rowSums(exp(temp))
  
  ##########################################
  # Remove negligible weights
  ##########################################
  
  keep <- w >=
    10 * .Machine$double.eps
  
  w <- w[keep]
  
  bs <- bs[keep]
  
  w <- w / sum(w)
  
  ##########################################
  # Posterior mean b
  ##########################################
  
  b <- sum(bs * w)
  
  ##########################################
  # Estimate k
  ##########################################
  
  temp <- log1p(-b * x)
  
  k <- mean(temp)
  
  ##########################################
  # Estimate sigma
  ##########################################
  
  sigma <- -k / b
  
  ##########################################
  # Weak prior shrinkage
  ##########################################
  
  a <- 10
  
  k <- (n * k + 0.5 * a) /
    (n + a)
  
  list(
    k = k,
    sigma = sigma
  )
}


gpinv <- function(p,
                  k,
                  sigma){
  
  ##########################################
  # Initialize output
  ##########################################
  
  x <- rep(NaN, length(p))
  
  ##########################################
  # Invalid sigma
  ##########################################
  
  if(sigma <= 0)
    return(x)
  
  ##########################################
  # Valid probabilities
  ##########################################
  
  ok <- (p > 0) & (p < 1)
  
  ##########################################
  # All probabilities valid
  ##########################################
  
  if(all(ok)){
    
    ########################################
    # k ≈ 0 (Exponential limit)
    ########################################
    
    if(abs(k) <
       .Machine$double.eps){
      
      x <- -log1p(-p)
      
    } else {
      
      ######################################
      # General GPD inverse
      ######################################
      
      x <-
        expm1(
          -k * log1p(-p)
        ) / k
    }
    
    x <- sigma * x
    
  } else {
    
    ########################################
    # Only subset valid
    ########################################
    
    p_ok <- p[ok]
    
    if(abs(k) <
       .Machine$double.eps){
      
      temp <- -log1p(-p_ok)
      
    } else {
      
      temp <-
        expm1(
          -k * log1p(-p_ok)
        ) / k
    }
    
    x[ok] <- sigma * temp
    
    ########################################
    # Boundary values
    ########################################
    
    x[p == 0] <- 0
    
    if(k >= 0){
      
      x[p == 1] <- Inf
      
    } else {
      
      x[p == 1] <- -sigma / k
    }
  }
  
  return(x)
}


psis_weights <- function(log_ratios,
                         Reff = 1,
                         tail_frac = 0.2){
  
  if(any(!is.finite(log_ratios)))
    stop("Non-finite log ratios detected.")
  
  S <- length(log_ratios)
  
  if(S <= 1)
    stop("More than one log-weight needed.")
  
  # Improve numerical accuracy
  
  lw <- log_ratios - max(log_ratios)
  
  # Random thresholding according to Vehtari etal
  
  M <- -ceiling(min(tail_frac * S, 3 * sqrt(S / Reff))) - 1
  
  cutoffmin <- log(.Machine$double.xmin)
  
  ord <- order(lw)
  
  xcutoff <- max(lw[ord][S + M + 1], cutoffmin)
  
  expxcutoff <- exp(xcutoff)
  

  # Tail samples
  
  tailinds <- which(lw > xcutoff)
  
  x2 <- lw[tailinds]
  
  n2 <- length(x2)
  
  # Not enough tail points
  
  if(n2 <= 4){
    
    k_hat <- Inf
    
  } else {
    
    # Sort tail
    
    x2si <- order(x2)
    
    # Excesses above threshold
    
    x2_exp <- exp(x2)
    
    x2_excess <- x2_exp - expxcutoff
    
    
    # Zhang–Stephens GPD fit
    
    fit <- gpdfitnew(x2_excess, sort = x2si)
    
    k_hat <- fit$k
    sigma_hat <- fit$sigma
  }
  
  # Pareto smoothing
  
  k_min <- 1 / 3
  
  if(is.finite(k_hat) && k_hat >= k_min){
    
    sti <- ((1:n2) - 0.5) / n2
    
    qq <- gpinv(sti, k_hat, sigma_hat)
    
    qq <- log(qq + expxcutoff)
    
    # Replace tail
    
    lw[tailinds[x2si]] <- qq
    
    # Truncate to largest raw weight = 0
    
    lw[lw > 0] <- 0
  }
  
  # Normalize in log-space
  
  lw <- lw - log_sum_exp(lw)
 
  w <- exp(lw)
  
  attr(w, "pareto_k") <- k_hat
  
  return(w)
}

tuple_function <- function(
    log_lik_mat,
    method,
    weight_type = "is",
    alpha = NULL
){
  
  method <- match.arg(
    method,
    c(
      "Gelfand",
      "Bronze",
      "Zanella"
    )
  )
  
  weight_type <- match.arg(
    weight_type,
    c(
      "is",
      "tis",
      "psis"
    )
  )
  
  S <- nrow(log_lik_mat)
  
  n <- ncol(log_lik_mat)
  
  ##########################################
  # Construct raw weights
  ##########################################
  
  w_mat <- matrix(
    0,
    S,
    n
  )
  
  ##########################################
  # Gelfand
  ##########################################
  
  if(method == "Gelfand"){
    
    w_mat <-
      exp(-log_lik_mat)
  }
  
  ##########################################
  # Bronze
  ##########################################
  
  if(method == "Bronze"){
    
    if(is.null(alpha)){
      stop(
        "alpha must be supplied for Bronze."
      )
    }
    
    w_mat <-
      exp(
        -alpha *
          log_lik_mat
      )
  }
  
  ##########################################
  # Zanella
  ##########################################
  
  if(method == "Zanella"){
    
    log_inv <- -log_lik_mat
    
    log_norm <- apply(
      log_inv,
      1,
      log_sum_exp
    )
    
    for(i in 1:n){
      
      w_mat[, i] <-
        exp(
          log_inv[, i] -
            log_norm
        )
    }
  }
  
  ##########################################
  # Weight stabilization
  ##########################################
  
  for(i in 1:n){
    
    ########################################
    # Ordinary IS
    ########################################
    
    if(weight_type == "is"){
      
      w <- w_mat[, i]
    }
    
    ########################################
    # Truncated IS
    ########################################
    
    if(weight_type == "tis"){
      
      w <-
        truncate_weights(
          w_mat[, i]
        )
    }
    
    ########################################
    # PSIS
    ########################################
    
    if(weight_type == "psis"){
      
      w <-
        psis_weights(
          log(w_mat[, i])
        )
    }
    
    w_mat[, i] <- w
  }
  
  ##########################################
  # Uniform tuple construction
  ##########################################
  
  numerator <-
    exp(log_lik_mat) *
    w_mat
  
  denominator <- w_mat
  
  list(
    numerator = numerator,
    denominator = denominator
  )
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
  
  if(is.null(num)){
    
    v <- -1 / colMeans(den)
    Y <- den
    mu <- 1 / colMeans(den)
    
  } else {
    
    v <- c(1 / colMeans(num), -1 / colMeans(den))
    Y <- cbind(num, den)
    mu <- colMeans(num) / colMeans(den)
    
  }
  
  qf <- bm.quad.form(Y = Y, v = v, b = b)
  
  list(elppd = sum(log(mu)), asymp_var = qf / S)
}


# psis_weights <- function(log_ratios,
#                          tail_frac = 0.2) {
#   
#   ############################################
#   # Remove non-finite values
#   ############################################
#   
#   keep <- is.finite(log_ratios)
#   
#   if(!all(keep))
#     stop("Non-finite log ratios detected.")
#   
#   ############################################
#   # Number of samples
#   ############################################
#   
#   S <- length(log_ratios)
#   
#   if(S < 10)
#     stop("Need at least 10 samples.")
#   
#   ############################################
#   # Stabilize
#   ############################################
#   
#   lw <- log_ratios - max(log_ratios)
#   
#   ############################################
#   # Raw importance ratios
#   ############################################
#   
#   rs <- exp(lw)
#   
#   ############################################
#   # Sort ratios
#   ############################################
#   
#   ord <- order(rs)
#   
#   rs_sorted <- rs[ord]
#   
#   ############################################
#   # Tail size
#   ############################################
#   
#   M <- floor(min(tail_frac * S,
#                  3 * sqrt(S)))
#   
#   M <- max(M, 5)
#   
#   ############################################
#   # Tail threshold
#   ############################################
#   
#   threshold_index <- S - M + 1
#   
#   u_hat <- rs_sorted[threshold_index]
#   
#   ############################################
#   # Excesses above threshold
#   ############################################
#   
#   excess <- rs_sorted[threshold_index:S] - u_hat
#   
#   ############################################
#   # Method-of-moments GPD estimates
#   ############################################
#   
#   m1 <- mean(excess)
#   
#   m2 <- mean(excess^2)
#   
#   denom <- m2 - 2 * m1^2
#   
#   if(denom <= 0 || !is.finite(denom)) {
#     
#     k_hat <- 0
#     
#     sigma_hat <- max(m1, 1e-8)
#     
#   } else {
#     
#     k_hat <- 0.5 * (m2 / denom - 1)
#     
#     sigma_hat <- 0.5 * m1 * (m2 / denom + 1)
#     
#   }
#   
#   ############################################
#   # Safeguards
#   ############################################
#   
#   if(!is.finite(k_hat))
#     k_hat <- 0
#   
#   if(!is.finite(sigma_hat) ||
#      sigma_hat <= 0)
#     sigma_hat <- max(m1, 1e-8)
#   
#   ############################################
#   # Pareto smoothing
#   ############################################
#   
#   ws_sorted <- rs_sorted
#   
#   p_grid <- ((1:M) - 0.5) / M
#   
#   for(z in 1:M) {
#     
#     p <- p_grid[z]
#     
#     ##########################################
#     # Smoothed quantile
#     ##########################################
#     
#     if(abs(k_hat) < 1e-8) {
#       
#       smoothed <- u_hat -
#         sigma_hat * log(1 - p)
#       
#     } else {
#       
#       smoothed <- u_hat +
#         (sigma_hat / k_hat) *
#         ((1 - p)^(-k_hat) - 1)
#     }
#     
#     ##########################################
#     # Replace tail
#     ##########################################
#     
#     ws_sorted[threshold_index + z - 1] <- min(
#       smoothed,
#       max(rs)
#     )
#   }
#   
#   ############################################
#   # Restore original order
#   ############################################
#   
#   ws <- numeric(S)
#   
#   ws[ord] <- ws_sorted
#   
#   ############################################
#   # Truncation
#   ############################################
#   
#   tau <- sqrt(S) * mean(ws)
#   
#   ws <- pmin(ws, tau)
#   
#   ############################################
#   # Normalize
#   ############################################
#   
#   ws <- ws / sum(ws)
#   
#   ############################################
#   # Pareto diagnostic warning
#   ############################################
#   
#   threshold <- min(1 - 1/log10(S),
#                    0.7)
#   
#   if(S > 2000)
#     threshold <- 0.7
#   
#   if(k_hat > threshold) {
#     
#     warning(
#       paste(
#         "Estimated Pareto k =",
#         round(k_hat, 3),
#         "suggests unstable importance sampling."
#       )
#     )
#   }
#   
#   return(ws)
# }
# 
# 
# # ELPPD tuple construction
# 
# tuple_function <- function(log_lik_mat,
#                            method,
#                            weight_type = "is",
#                            alpha = NULL){
#   
#   method <- match.arg(method,
#                       c("Gelfand",
#                         "Bronze",
#                         "Zanella"))
#   
#   weight_type <- match.arg(weight_type,
#                            c("is",
#                              "tis",
#                              "psis"))
#   
#   S <- nrow(log_lik_mat)
#   
#   n <- ncol(log_lik_mat)
#   
#   
#   # ---------------------------------------------------------
#   # Classical importance weights
#   # ---------------------------------------------------------
#   
#   if(method == "Gelfand"){
#     
#     numerator <- matrix(1, S, n)
#     
#     denominator <- exp(-log_lik_mat)
#   }
#   
#   
#   if(method == "Bronze"){
#     
#     numerator <- exp((1 - alpha) * log_lik_mat)
#     
#     denominator <- exp(-alpha * log_lik_mat)
#   }
#   
#   
#   if(method == "Zanella"){
#     
#     log_inv <- -log_lik_mat
#     
#     log_denom <- apply(log_inv, 1, log_sum_exp)
#     
#     denominator <- matrix(0, S, n)
#     
#     for(i in 1:n){
#       
#       denominator[, i] <- exp(log_inv[, i] - log_denom)
#     }
#     
#     numerator <- exp(log_lik_mat) * denominator
#   }
#   
#   
#   # ---------------------------------------------------------
#   # Weight stabilization
#   # ---------------------------------------------------------
#   
#   if(weight_type == "tis"){
#     
#     for(i in 1:n){
#       
#       denominator[, i] <- truncate_weights(denominator[, i])
#     }
#   }
#   
#   
#   if(weight_type == "psis"){
#     
#     for(i in 1:n){
#       
#       denominator[, i] <- psis_weights(
#         log(denominator[, i])
#       )
#     }
#   }
#   
#   
#   list(
#     numerator = numerator,
#     denominator = denominator
#   )
# }
