tuple_function <- function(log_lik_mat,
                           method,
                           alpha_obs = NULL){
  
  method <- match.arg(method, c("Gelfand", "Bronze", "Zanella"))
  
  S <- nrow(log_lik_mat)
  
  n <- ncol(log_lik_mat)
  
  # Gelfand
  
  if(method == "Gelfand"){
    
    numerator <- matrix(
      1,
      nrow(log_lik_mat),
      ncol(log_lik_mat)
    )
    
    denominator <- exp(-log_lik_mat)
  }
  

  # Bronze
  
  if(method == "Bronze"){
    
    numerator <- exp(log_lik_mat * matrix(1 - alpha_obs, 
                                          nrow(log_lik_mat),
                                          ncol(log_lik_mat),
                                          byrow = TRUE
                                          ))
    
    denominator <- exp(log_lik_mat * matrix(-alpha_obs,
                                            nrow(log_lik_mat),
                                            ncol(log_lik_mat),
                                            byrow = TRUE))
  }
  
  
  # Zanella
  
  if(method == "Zanella"){
    
    log_inv <- -log_lik_mat
    
    log_denom <- apply(
      log_inv,
      1,
      log_sum_exp
    )
    
    denominator <- matrix(
      0,
      S,
      n
    )
    
    for(i in 1:n){
      
      denominator[, i] <- exp(
        log_inv[, i] -
          log_denom
      )
    }
    
    numerator <- exp(log_lik_mat) *
      denominator
  }
  
  list(
    numerator = numerator,
    denominator = denominator
  )
}
# Truncated IS


truncate_weights <- function(w){
  
  tau <- sqrt(length(w)) * mean(w)
  
  pmin(w, tau)
}

bm.chol.omega <- function(Y,
                          b = NULL){
  
  Y <- as.matrix(Y)
  
  S <- nrow(Y)
  
  p <- ncol(Y)
  
  if(is.null(b)){
    
    b <- floor(sqrt(S))
  }
  
  a <- floor(S / b)
  
  S_trim <- a * b
  
  Y <- Y[
    1:S_trim,
    ,
    drop = FALSE
  ]
  
  batch_means <- matrix(0, a, p)
  
  for(k in 1:a){
    
    idx <- ((k - 1) * b + 1):(k * b)
    
    batch_means[k, ] <-
      colMeans(
        Y[idx, , drop = FALSE]
      )
  }
  
  grand_mean <- colMeans(batch_means)
  
  B <- sqrt(b / (a - 1)) *
    sweep(
      batch_means,
      2,
      grand_mean
    )
  
  chol(crossprod(B))
}

bm.quad.form <- function(Y,
                         v,
                         b = NULL){
  
  U <- bm.chol.omega(Y, b)
  
  Uv <- U %*% v
  
  as.numeric(crossprod(Uv))
}

asymp.var.elppd <- function(num,
                            den,
                            S,
                            b = NULL){
  
  if(is.null(num)){
    
    v <- -1 / colMeans(den)
    
    Y <- den
    
    mu <- 1 / colMeans(den)
    
  } else {
    
    v <- c(
      1 / colMeans(num),
      -1 / colMeans(den)
    )
    
    Y <- cbind(num, den)
    
    mu <- colMeans(num) /
      colMeans(den)
  }
  
  qf <- bm.quad.form(
    Y = Y,
    v = v,
    b = b
  )
  
  list(
    elppd = sum(log(mu)),
    asymp_var = qf / S
  )
}


