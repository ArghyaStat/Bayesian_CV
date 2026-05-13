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

gpd_quantile <- function(p, k, sigma, u) {
  
  if(abs(k) < 1e-8) {
    
    return(
      u - sigma * log(1 - p)
    )
    
  } else {
    
    return(
      u + (sigma / k) * ((1 - p)^(-k) - 1)
    )
    
  }
  
}


psis_weights <- function(log_ratios) {
  
  S <- length(log_ratios)
  
  ############################################
  # Stabilize
  ############################################
  
  log_ratios <- log_ratios - max(log_ratios)
  
  ############################################
  # Raw importance ratios
  ############################################
  
  rs <- exp(log_ratios)
  
  ############################################
  # Order statistics
  ############################################
  
  ord <- order(rs)
  
  rs_sorted <- rs[ord]
  
  ############################################
  # Tail size
  ############################################
  
  M <- floor(min(0.2 * S, 3 * sqrt(S)))
  
  M <- max(M, 5)
  
  ############################################
  # Threshold:
  # T = r_{S-M+1:S}
  ############################################
  
  threshold_index <- S - M + 1
  
  u_hat <- rs_sorted[threshold_index]
  
  ############################################
  # Tail excesses
  ############################################
  
  excess <- rs_sorted[threshold_index:S] - u_hat
  
  ############################################
  # Method-of-moments GPD estimates
  ############################################
  
  m1 <- mean(excess)
  
  m2 <- mean(excess^2)
  
  denom <- m2 - m1^2
  
  if(denom <= 0 || !is.finite(denom)) {
    
    k_hat <- 0
    
    sigma_hat <- m1
    
  } else {
    
    k_hat <- 0.5 * (m2 / denom - 1)
    
    sigma_hat <- 0.5 * m1 * (m2 / denom + 1)
    
  }
  
  ############################################
  # Numerical safeguards
  ############################################
  
  if(!is.finite(k_hat))
    k_hat <- 0
  
  if(!is.finite(sigma_hat) || sigma_hat <= 0)
    sigma_hat <- m1
  
  ############################################
  # Pareto smoothing of upper tail
  ############################################
  
  ws_sorted <- rs_sorted
  
  for(z in 1:M) {
    
    p <- (z - 0.5) / M
    
    smoothed <- gpd_quantile(
      p = p,
      k = k_hat,
      sigma = sigma_hat,
      u = u_hat
    )
    
    ##########################################
    # Truncate at maximum raw ratio
    ##########################################
    
    ws_sorted[threshold_index + z - 1] <- min(
      smoothed,
      max(rs)
    )
    
  }
  
  ############################################
  # Restore original order
  ############################################
  
  ws <- numeric(S)
  
  ws[ord] <- ws_sorted
  
  ############################################
  # Self-normalized IS weights
  ############################################
  
  ws <- ws / sum(ws)
  
  ############################################
  # Pareto-k diagnostic
  ############################################
  
  threshold <- min(1 - 1/log10(S), 0.7)
  
  if(S > 2000)
    threshold <- 0.7
  
  if(k_hat > threshold) {
    
    warning(
      paste(
        "Estimated Pareto k =",
        round(k_hat, 3),
        "suggests unstable importance sampling."
      )
    )
    
  }
  
  return(ws)
  
}



# ELPPD tuple construction

tuple_function <- function(log_lik_mat,
                           method,
                           weight_type = "is",
                           alpha = NULL){
  
  method <- match.arg(method,
                      c("Gelfand",
                        "Bronze",
                        "Zanella"))
  
  weight_type <- match.arg(weight_type,
                           c("is",
                             "tis",
                             "psis"))
  
  S <- nrow(log_lik_mat)
  
  n <- ncol(log_lik_mat)
  
  
  # ---------------------------------------------------------
  # Classical importance weights
  # ---------------------------------------------------------
  
  if(method == "Gelfand"){
    
    numerator <- matrix(1, S, n)
    
    denominator <- exp(-log_lik_mat)
  }
  
  
  if(method == "Bronze"){
    
    numerator <- exp((1 - alpha) * log_lik_mat)
    
    denominator <- exp(-alpha * log_lik_mat)
  }
  
  
  if(method == "Zanella"){
    
    log_inv <- -log_lik_mat
    
    log_denom <- apply(log_inv, 1, log_sum_exp)
    
    denominator <- matrix(0, S, n)
    
    for(i in 1:n){
      
      denominator[, i] <- exp(log_inv[, i] - log_denom)
    }
    
    numerator <- exp(log_lik_mat) * denominator
  }
  
  
  # ---------------------------------------------------------
  # Weight stabilization
  # ---------------------------------------------------------
  
  if(weight_type == "tis"){
    
    for(i in 1:n){
      
      denominator[, i] <- truncate_weights(denominator[, i])
    }
  }
  
  
  if(weight_type == "psis"){
    
    for(i in 1:n){
      
      denominator[, i] <- psis_weights(
        denominator[, i]
      )
    }
  }
  
  
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