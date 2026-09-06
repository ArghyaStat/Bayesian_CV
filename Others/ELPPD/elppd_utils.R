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
  
  method <- match.arg(method, c("Gelfand", "Bronze", "Zanella"))
  
  weight_type <- match.arg(weight_type, c("is", "tis", "psis"))
  
  S <- nrow(log_lik_mat)
  n <- ncol(log_lik_mat)
  
  if(length(kappa) != n){
    stop("length(kappa) must equal ncol(log_lik_mat)")
  }
  
  K <- max(kappa)
  
  # Fold log-likelihoods
  
  fold_loglik <- matrix(0, S, K)
  
  for(k in seq_len(K)){
    
    idx <- which(kappa == k)
    
    # kth column chooses those indices idx which are in k-th fold.
    
    fold_loglik[, k] <- rowSums(log_lik_mat[, idx, drop = FALSE])
  }
  
  # total_loglik is a S length vector 
  # computing total likelihood from n observations
  
  total_loglik <- rowSums(log_lik_mat)
  
  # Weight matrix
  
  w_mat <- matrix(0, S, n)
  
  for(k in seq_len(K)){
    
    idx <- which(kappa == k)
    
    # Gelfand
    
    if(method == "Gelfand"){
      
      log_w_k <- -fold_loglik[, k]
    }
    
    # Bronze
    
    if(method == "Bronze"){
      
      if(is.null(alpha)){
        stop("alpha must be supplied for Bronze.")
      }
      
      log_w_k <- (1 - alpha) * total_loglik - fold_loglik[, k]
    }
    
    # Zanella
    
    if(method == "Zanella"){
      
      # 1/prod_k(p(y_j | \theta))
      
      log_inv_fold <- -fold_loglik
      
      log_norm <- apply(log_inv_fold, 1, log_sum_exp)
      
      log_w_k <- log_inv_fold[, k] - log_norm
    }
    
    # numerical stabilization only
    log_w_k <- log_w_k - max(log_w_k)
    
    # IS / TIS / PSIS
    if(weight_type == "is"){
      
      w_k <- exp(log_w_k)
    }
    
    if(weight_type == "tis"){
      
      w_k <- truncate_weights(exp(log_w_k))
      
    }
    
    if(weight_type == "psis"){
      
      r_eff_k <- loo::relative_eff(exp(log_w_k), chain_id = rep(1, S))
      psis_details <- loo::psis(log_w_k, r_eff = r_eff_k)

      w_k <- weights(psis_details, log = FALSE, normalize = FALSE)
    }
    
    w_mat[, idx] <- matrix(w_k, nrow = S, ncol = length(idx))
  
    }
    
  
  # Tuple construction
  
  numerator <- exp(log_lik_mat) * w_mat
  
  denominator <- w_mat
  
  list(numerator = numerator,
       denominator = denominator)
}


# Ordinary batch means estimator

bm.var.base <- function(x, b){
  
  x <- as.numeric(x)
  S <- length(x)
  
  if(b < 2){
    stop("Batch size b must be >= 2.")
  }
  
  a <- floor(S / b)
  
  if(a < 2){
    stop("Need at least two batches.")
  }
  
  S_trim <- a * b
  
  x <- x[seq_len(S_trim)]
  
  dim(x) <- c(b, a)
  
  batch_means <- colMeans(x)
  
  grand_mean <- mean(batch_means)
  
  sigma2_hat <- b * sum((batch_means - grand_mean)^2) /(a - 1)
  
  sigma2_hat
}


# Lugsail batch means estimator


bm.var.ls <- function(x, b = NULL, r = 3){
  
  x <- as.numeric(x)
  
  S <- length(x)
  
  if(S < 10){
    stop("Chain too short.")
  }
  
  # default batch size
  if(is.null(b)){
    b <- floor(sqrt(S))
  }
  
  if(r < 1){
    stop("r must be >= 1.")
  }
  
  # ordinary BM (r = 1)
  if(r == 1){
    return(bm.var.base(x, b))
  }
  
  # main batch size
  sigma2_b <- bm.var.base(x, b)
  
  # reduced batch size
  b_small <- max(2, floor(b / r))
  
  sigma2_small <- bm.var.base(x, b_small)
  
  # lugsail correction
  sigma2_lugsail <-(r / (r - 1)) * sigma2_b - (1 / (r - 1)) * sigma2_small
  
  sigma2_lugsail
}

# Asymptotic variance estimator

elppd.details <- function(num, 
                          den, 
                          S = NULL, 
                          b = NULL,
                          method = c("our", "mcmcse")){
  
  method <- match.arg(method)
  num <- as.matrix(num)
  den <- as.matrix(den)
  
  mu_num <- colMeans(num)
  mu_den <- colMeans(den)
  
  Y <- cbind(num, den)
  S <- nrow(Y)
  n <- ncol(num)
  
  se2_num <- numeric(n)
  se2_den <- numeric(n)
  
  if (method == "our") {
    for (i in seq_len(n)) {
      se2_num[i] <- bm.var.ls(num[, i], b = b, r = 3) / S
      se2_den[i] <- bm.var.ls(den[, i], b = b, r = 3) / S
    }
  } else if (method == "mcmcse") {
    for (i in seq_len(n)) {
      se2_num[i] <- mcmcse::mcse(num[, i])$se^2
      se2_den[i] <- mcmcse::mcse(den[, i])$se^2
    }
  }
  
  v_num <- mu_num / (mu_num^2 + se2_num)
  v_den <- -mu_den / (mu_den^2 + se2_den)
  v <- c(v_num, v_den)
  
  if(S < 2){
    stop("Need at least two MCMC draws.")
  }
  
  # ELPPD estimate
  
  mu <- mu_num / mu_den
  elppd_hat <- sum(log(mu))
  
  xi <- as.numeric(Y %*% v)
  
  
  # Note: The quantities below estimate the variance of the estimator, i.e.,
  # Var(ELPPD_hat) = sigma2_ELPPD / S, where sigma2_ELPPD is the asymptotic
  # variance from the CLT, and S is the number of MCMC draws.
  if(method == "our"){
    
    asymp_se <- bm.var.ls(xi, b = NULL, r = 3)/S
  }
  
  if(method == "mcmcse"){
    
    asymp_se <- (mcmcse::mcse(xi)$se^2)

    }
  
  list(
    elppd = elppd_hat,
    asymp_se =  asymp_se)
}
