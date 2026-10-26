# Vecchia ordering and conditioning sets

# Max-min ordering (Guinness, 2018)
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

# m nearest previously-ordered neighbours of each location (row includes the location itself)
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
