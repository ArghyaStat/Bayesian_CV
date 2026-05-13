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

# Generalized Pareto quantile function

qgpd <- function(p, k, sigma){
  
  if(abs(k) < 1e-8){
    
    -sigma * log(1 - p)
    
  } else {
    
    sigma / k * ((1 - p)^(-k) - 1)
  }
}

# Pareto smoothed importance sampling weights

psis_weights <- function(w,
                         tail_frac = 0.2){
  
  S <- length(w)
  
  # Number of tail samples
  
  M <- ceiling(tail_frac * S)
  
  # Order statistics
  
  ord <- order(w)
  
  w_sorted <- w[ord]
  
  # Threshold
  
  threshold_idx <- S - M + 1
  
  threshold <- w_sorted[threshold_idx]
  
  # Tail excesses
  
  tail_weights <- w_sorted[threshold_idx:S]
  
  excesses <- tail_weights - threshold
  
  # Fit generalized Pareto to excesses
  
  fit <- loo::gpdfit(excesses)
  
  k <- fit$k
  sigma <- fit$sigma
  
  # Smoothed tail quantiles
  
  p <- ((1:M) - 0.5) / M
  
  smoothed_excess <- qgpd(p,
                          k = k,
                          sigma = sigma)
  
  smoothed_tail <- threshold + smoothed_excess
  
  # Replace tail weights
  
  w_smoothed <- w_sorted
  
  w_smoothed[threshold_idx:S] <- smoothed_tail
  
  # Optional stabilization:
  # truncate at largest raw weight
  
  w_smoothed <- pmin(w_smoothed,
                     max(w))
  
  # Return original ordering
  
  out <- numeric(S)
  
  out[ord] <- w_smoothed
  
  out
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
      
      denominator[, i] <- psis_weights(denominator[, i])
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