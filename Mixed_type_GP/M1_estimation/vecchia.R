
### Fast maxmin ordering (Guiness 2018)

max_min <- function(locations) {
  n <- nrow(locations)
  m <- round(sqrt(n))
  
  NNall <- get.knn(locations, k = m)$nn.index
  
  centroid <- colMeans(locations)
  dist_centroid <- sqrt(rowSums((locations - matrix(centroid,n,ncol(locations),TRUE))^2))
  first_point <- which.min(dist_centroid)
  
  ordered_indices <- integer(n)
  ordered_indices[1] <- first_point
  
  selected <- rep(FALSE,n)
  selected[first_point] <- TRUE
  
  min_dists <- rep(Inf,n)
  
  for(i in seq_len(n)){
    min_dists[i] <- sqrt(sum((locations[i,]-locations[first_point,])^2))
  }
  
  for(k in 2:n){
    
    next_point <- which.max(min_dists)
    
    ordered_indices[k] <- next_point
    selected[next_point] <- TRUE
    
    min_dists[next_point] <- -Inf
    
    for(i in which(!selected)){
      
      new_dist <- sqrt(sum((locations[i,]-locations[next_point,])^2))
      min_dists[i] <- min(min_dists[i], new_dist)
      
    }
    
  }
  
  ordered_indices
}

### preceding ordered neighbouring

neighbor_matrix <- function(locs_ord, m) {
  
  n <- nrow(locs_ord)
  
  nn_matrix <- matrix(NA, nrow = n, ncol = m + 1)
  
  for (i in 1:n) {
    
    pred_ind <- 1:i
    
    dx <- locs_ord[pred_ind, 1] - locs_ord[i, 1]
    dy <- locs_ord[pred_ind, 2] - locs_ord[i, 2]
    
    dists <- dx*dx + dy*dy   # squared distance (faster than sqrt)
    
    nearest_indices_within_subset <- order(dists)
    
    sorted_pred_ind <- pred_ind[nearest_indices_within_subset]
    
    sorted_neighbors <- sorted_pred_ind[1:min(m + 1, length(sorted_pred_ind))]
    
    row <- rep(NA, m + 1)
    row[1:length(sorted_neighbors)] <- sort(sorted_neighbors)
    
    nn_matrix[i, ] <- row
  }
  
  row.names(nn_matrix) <- NULL
  
  return(nn_matrix)
}



# Distance computation function returning a list of n matrices
dist.nn <- function(locs.ord, neighbor_matrix) {
  
  n <- nrow(locs.ord)
  
  dist_nn <- list()
  for (i in 1:n) {
    neighbors <- neighbor_matrix[i, ]
    neighbors <- neighbors[!is.na(neighbors)]  # Remove NA values
    
    if (length(neighbors) > 0) {
      dist_nn[[i]] <- rdist(locs.ord[neighbors, , drop = FALSE], locs.ord[neighbors, , drop = FALSE])
    } else {
      dist_nn[[i]] <- NULL
    }
  }
  
  return(dist_nn)
}


### Vecchia_coefficcients for faster computation

vecchia_coeff <- function(dist.nn, NNarray, phi, nu){
  
  n <- length(dist.nn)
  
  B_list <- vector("list", n)
  r_vec  <- numeric(n)
  
  for(i in seq_len(n)){
    
    neigh <- NNarray[i, ]
    neigh <- neigh[!is.na(neigh)]
    
    d <- dist.nn[[i]]
    
    K <- Matern(d, range = phi, smoothness = nu)
    
    # find where location i appears inside the block
    idx_i <- which(neigh == i)
    
    if(length(neigh) == 1){
      
      r_vec[i] <- K[idx_i, idx_i]
      B_list[[i]] <- numeric(0)
      
    } else {
      
      K11 <- K[idx_i, idx_i]
      
      K12 <- K[idx_i, -idx_i]
      K21 <- K[-idx_i, idx_i]
      K22 <- K[-idx_i, -idx_i, drop = FALSE]
      
      R <- chol(K22)
      
      x <- backsolve(R, forwardsolve(t(R), K21))
      
      B_i <- as.numeric(t(x))
      
      r_vec[i] <- K11 - sum(B_i * K12)
      
      B_list[[i]] <- B_i
    }
  }
  
  list(B = B_list, r = r_vec)
}



### Fast vecchia likelihood


vecchia_ll <- function(W.ord, M.ord, NNarray, coeff, chol.V, log = TRUE){
  
  B_list <- coeff$B
  r_vec  <- coeff$r
  
  n <- nrow(W.ord)
  q <- ncol(W.ord)
  
  # log|Sigma|
  logdetV <- 2 * sum(log(diag(chol.V)))
  
  ll <- 0
  
  for(i in seq_len(n)){
    
    Xi <- W.ord[i,,drop=FALSE]
    Mi <- M.ord[i,,drop=FALSE]
    
    neigh <- NNarray[i,]
    neigh <- neigh[!is.na(neigh)]
    
    idx_i <- which(neigh == i)
    neigh_cond <- neigh[-idx_i]
    
    if(length(neigh_cond) == 0){
      
      resid <- Xi - Mi
      
    } else {
      
      Xn <- W.ord[neigh_cond,,drop=FALSE]
      Mn <- M.ord[neigh_cond,,drop=FALSE]
      
      B <- matrix(B_list[[i]],1)
      
      resid <- Xi - Mi - B %*% (Xn - Mn)
    }
    
    r <- r_vec[i]
    
    resid <- as.numeric(resid)
    
    z <- forwardsolve(t(chol.V), resid)
    
    quad <- sum(z^2) / r
    
    # u <- resid %*% t(chol.V.prec)
    # quad <- sum(u^2) / r
    
    ll <- ll -0.5 * ( q*log(2*pi*r) + logdetV + quad )
  }
  
  if(log) return(as.numeric(ll)) else return(exp(ll))
}

## Vecchia products UW and UX

vecchia_products <- function(W, X, B, r, NNarray){
  
  N <- nrow(W)
  q <- ncol(W)
  p <- ncol(X)
  
  UW <- matrix(0, N, q)
  UX <- matrix(0, N, p)
  
  
  for(i in seq_len(N)){
    
    neigh <- NNarray[i,]
    neigh <- neigh[!is.na(neigh)]
    neigh <- neigh[neigh < i]   
    
    wi <- W[i,]
    xi <- X[i,]
    
    if(length(neigh) > 0){
      
      B_i <- B[[i]]
      
      ## ---- CRITICAL FIX: match lengths ----
      B_i <- B_i[seq_along(neigh)]
      
      ## weighted sums (no matrix mult → no conformability issues)
      wi <- wi - colSums(W[neigh,,drop=FALSE] * B_i)
      xi <- xi - colSums(X[neigh,,drop=FALSE] * B_i)
    }
    
   
    
    UW[i,] <- wi/sqrt(r[i])
    UX[i,] <- xi/sqrt(r[i])
  }
  
  list(
    UW = UW,
    UX = UX
  )
}




vecchia_prior_draw <- function(Z, NNarray, B_list, r_vec){

  N <- nrow(Z)
  q <- ncol(Z)

  nu <- matrix(0, N, q)

  for(i in seq_len(N)){

    neigh <- NNarray[i,]
    neigh <- neigh[!is.na(neigh)]
    neigh <- neigh[neigh < i]

    if(length(neigh) == 0){

      nu[i,] <- sqrt(r_vec[i]) * Z[i,]

    } else {

      B <- B_list[[i]]
      nu[i,] <- as.numeric(B %*% nu[neigh,,drop=FALSE]) +
        sqrt(r_vec[i]) * Z[i,]

    }
  }

  nu
}


vecchia_pred_draw <- function(W.obs.ord,
                              beta,
                              Sigma,
                              X.obs.ord,
                              X.pred.ord,
                              NNarray,
                              coeff,
                              obs_id,
                              pred_id){
  
  B_list <- coeff$B
  r_vec  <- coeff$r
  
  q <- ncol(W.obs.ord)
  n_tot <- nrow(NNarray)
  
  mu <- matrix(0, n_tot, q)
  mu[obs_id,]  <- X.obs.ord %*% beta
  mu[pred_id,] <- X.pred.ord %*% beta
  
  W <- matrix(0, n_tot, q)
  W[obs_id,] <- W.obs.ord
  
  R_Sigma <- chol(Sigma)
  
  ## pre-generate all noise
  z.mat <- matrix(rnorm(length(pred_id)*q), length(pred_id), q)
  
  for(k in seq_along(pred_id)){
    
    i <- pred_id[k]
    
    neigh <- NNarray[i,]
    neigh <- neigh[!is.na(neigh)]
    
    B_i <- B_list[[i]]
    
    if(length(B_i) == 0){
      
      mean_i <- mu[i,]
      
    } else {
      
      idx <- neigh[neigh < i]
      
      diff <- W[idx,,drop=FALSE] - mu[idx,,drop=FALSE]
      
      mean_i <- mu[i,] + as.numeric(B_i %*% diff)
    }
    
    
    W[i,] <- mean_i + sqrt(r_vec[i]) * (z.mat[k, ] %*% R_Sigma)
  }
  
  W[pred_id,,drop=FALSE]
}
