
# Gelfand

tuple_g <- function(log_lik_mat){
  
  den <- exp(-log_lik_mat)
  
  list(
    numerator = NULL,
    denominator = den
  )
}



# ---------------------------------------------------------
# Bronze
# ---------------------------------------------------------

tuple_b <- function(log_lik_mat,
                    alpha_obs){
  
  num <- exp(
    sweep(log_lik_mat,
          2,
          1 - alpha_obs,
          "*")
  )
  
  den <- exp(
    sweep(log_lik_mat,
          2,
          -alpha_obs,
          "*")
  )
  
  list(
    numerator = num,
    denominator = den
  )
}



# ---------------------------------------------------------
# Zanella
# ---------------------------------------------------------

tuple_z <- function(log_lik_mat){
  
  S <- nrow(log_lik_mat)
  
  n <- ncol(log_lik_mat)
  
  log_inv <- -log_lik_mat
  
  log_denom <- apply(
    log_inv,
    1,
    log_sum_exp
  )
  
  den <- matrix(0, S, n)
  
  for(i in 1:n){
    
    den[, i] <- exp(
      log_inv[, i] -
        log_denom
    )
  }
  
  num <- exp(log_lik_mat) * den
  
  list(
    numerator = num,
    denominator = den
  )
}

# Truncated IS
# ---------------------------------------------------------

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


