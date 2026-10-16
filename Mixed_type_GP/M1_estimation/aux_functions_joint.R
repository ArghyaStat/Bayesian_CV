######### Posterior Computations #########

log.likelihood <- function(
    W,
    Y,
    family,
    family.par = NULL,
    link = NULL
) {
  
  # -------------------------------------------------
  # Checks
  # -------------------------------------------------
  
  if (!all(dim(W) == dim(Y)))
    stop("Dimensions of W and Y must match.")
  
  q <- ncol(W)
  
  if (length(family) != q)
    stop("Length of 'family' must equal number of responses.")
  
  if (is.null(family.par))
    family.par <- vector("list", q)
  
  if (length(family.par) != q)
    stop("family.par must be list of length q.")
  
  if (is.null(link))
    link <- rep("canonical", q)
  
  if (length(link) != q)
    stop("Length of 'link' must equal q.")
  
  
  # -------------------------------------------------
  # Link inverse
  # -------------------------------------------------
  
  linkinv <- function(eta, fam, linkname) {
    
    if (linkname == "canonical") {
      
      switch(fam,
             "Gaussian" = eta,
             "Poisson"  = exp(eta),
             "Binomial" = plogis(eta),
             "Gamma"    = exp(eta),
             "Negative-Binomial" = exp(eta),
             stop(paste("Unsupported family:", fam)))
      
    } else if (linkname == "logit") {
      
      plogis(eta)
      
    } else if (linkname == "probit") {
      
      pnorm(eta)
      
    } else if (linkname == "cloglog") {
      
      1 - exp(-exp(eta))
      
    } else if (linkname == "identity") {
      
      eta
      
    } else if (linkname == "log") {
      
      exp(eta)
      
    } else {
      
      stop("Unsupported link.")
    }
  }
  
  
  # -------------------------------------------------
  # Compute log-likelihood
  # -------------------------------------------------
  
  like.component <- numeric(q)
  
  for (j in seq_len(q)) {
    
    fam  <- family[j]
    pars <- family.par[[j]]
    eta  <- W[, j]
    mu   <- linkinv(eta, fam, link[j])
    
    if (fam == "Gaussian") {
      
      dispersion <- ifelse(is.null(pars$dispersion), 1, pars$dispersion)
      
      like.component[j] <- sum(
        dnorm(Y[, j], mean = mu,
              sd = sqrt(dispersion),
              log = TRUE)
      )
      
    } else if (fam == "Poisson") {
      
      like.component[j] <- sum(
        dpois(Y[, j], lambda = mu, log = TRUE)
      )
      
    } else if (fam == "Binomial") {
      
      size <- ifelse(is.null(pars$size), 1, pars$size)
      
      like.component[j] <- sum(
        dbinom(Y[, j], size = size,
               prob = mu,
               log = TRUE)
      )
      
    } else if (fam == "Gamma") {
      
      shape <- ifelse(is.null(pars$shape), 1, pars$shape)
      rate  <- shape / mu
      
      like.component[j] <- sum(
        dgamma(Y[, j],
               shape = shape,
               rate = rate,
               log = TRUE)
      )
      
    } else if (fam == "Negative-Binomial") {
      
      theta <- ifelse(is.null(pars$theta), 1, pars$theta)
      
      like.component[j] <- sum(
        dnbinom(Y[, j],
                mu = mu,
                size = theta,
                log = TRUE)
      )
      
    } else {
      
      stop(paste("Unsupported family:", fam))
    }
  }
  
  list(
    like.total = sum(like.component),
    like.comp  = like.component
  )
}


Y.pred.ord <- function(W.pred.ord, N.pred, q, family, family.par = NULL, link = NULL) {
  
  # Checks
  
  if (!all(dim(W.pred.ord) == c(N.pred, q)))
    stop("W.pred.ord must be N.pred x q.")
  
  if (length(family) != q)
    stop("Length of 'family' must equal q.")
  
  if (is.null(family.par))
    family.par <- vector("list", q)
  
  if (length(family.par) != q)
    stop("family.par must be list of length q.")
  
  if (is.null(link))
    link <- rep("canonical", q)
  
  if (length(link) != q)
    stop("Length of 'link' must equal q.")
  
  # Link inverse
  
  linkinv <- function(eta, fam, linkname) {
    
    if (linkname == "canonical") {
      
      switch(fam,
             "Gaussian" = eta,
             "Poisson"  = exp(eta),
             "Binomial" = plogis(eta),
             "Gamma"    = exp(eta),
             "Negative-Binomial" = exp(eta),
             stop(paste("Unsupported family:", fam)))
      
    } else if (linkname == "logit") {
      
      plogis(eta)
      
    } else if (linkname == "probit") {
      
      pnorm(eta)
      
    } else if (linkname == "cloglog") {
      
      1 - exp(-exp(eta))
      
    } else if (linkname == "identity") {
      
      eta
      
    } else if (linkname == "log") {
      
      exp(eta)
      
    } else {
      
      stop("Unsupported link.")
    }
  }
  
  # Generate predictive responses
  
  Y.mat <- matrix(NA, N.pred, q)
  
  for (j in seq_len(q)) {
    
    fam  <- family[j]
    pars <- family.par[[j]]
    eta  <- W.pred.ord[, j]
    
    mu <- linkinv(eta, fam, link[j])
    
    if (fam == "Gaussian") {
      
      dispersion <- ifelse(is.null(pars$dispersion), 1, pars$dispersion)
      
      Y.mat[, j] <- rnorm(
        N.pred,
        mean = mu,
        sd   = sqrt(dispersion)
      )
      
    } else if (fam == "Poisson") {
      
      Y.mat[, j] <- rpois(
        N.pred,
        lambda = mu
      )
      
    } else if (fam == "Binomial") {
      
      size <- ifelse(is.null(pars$size), 1, pars$size)
      
      Y.mat[, j] <- rbinom(
        N.pred,
        size = size,
        prob = mu
      )
      
    } else if (fam == "Gamma") {
      
      shape <- ifelse(is.null(pars$shape), 1, pars$shape)
      rate  <- shape / mu
      
      Y.mat[, j] <- rgamma(
        N.pred,
        shape = shape,
        rate  = rate
      )
      
    } else if (fam == "Negative-Binomial") {
      
      theta <- ifelse(is.null(pars$theta), 1, pars$theta)
      
      Y.mat[, j] <- rnbinom(
        N.pred,
        mu   = mu,
        size = theta
      )
      
    } else {
      
      stop(paste("Unsupported family:", fam))
    }
  }
  
  Y.mat
}


# Stable log-mean-exp function
log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}



compute_ess <- function(chain) {
  # Check if the chain is a list
  if (is.list(chain)) {
    # If the chain is a list of matrices, compute ESS element-wise and keep matrix structure
    if (all(sapply(chain, is.matrix))) {
      n_iter <- length(chain)
      n_row <- nrow(chain[[1]])
      n_col <- ncol(chain[[1]])
      
      # Initialize a matrix to store ESS values for each element of the matrix chain
      ess_matrix <- matrix(NA, nrow = n_row, ncol = n_col)
      
      # Compute ESS for each element of the matrix
      for (i in 1:n_row) {
        for (j in 1:n_col) {
          # Extract the (i, j) element across all iterations and compute ESS
          element_chain <- sapply(chain, function(x) x[i, j])
          ess_matrix[i, j] <- ess(element_chain)
        }
      }
      
      return(ess_matrix)  # Return the ESS matrix
      
    } else if (all(sapply(chain, is.vector))) {
      # If the chain is a list of numeric vectors (e.g., phi.chain), compute ESS for each vector
      ess_vector <- sapply(chain, ess)  # Apply ESS to each vector in the list
      return(ess_vector)
      
    } else {
      stop("The chain contains elements that are neither matrices nor vectors.")
    }
    
  } else if (is.vector(chain)) {
    # If the chain is a single numeric vector (e.g., phi.chain), compute ESS for the entire vector
    return(ess(chain))  # Directly apply ESS to the vector
    
  } else {
    stop("Unexpected structure for chain")
  }
}


summary_stats <- function(chain, true_value) {
  
  n_iter <- length(chain)
  
  if (is.matrix(chain[[1]])) {
    
    d1 <- nrow(chain[[1]])
    d2 <- ncol(chain[[1]])
    
    #### stack chain ####
    arr <- simplify2array(chain)
    
    #### posterior mean ####
    post.mean <- apply(arr, c(1,2), mean)
    
    #### posterior median ####
    post.median <- apply(arr, c(1,2), median)
    
    #### full covariance of vec(parameter) ####
    theta <- matrix(NA, n_iter, d1*d2)
    
    for (i in seq_len(n_iter)) {
      theta[i,] <- as.vector(chain[[i]])
    }
    
    cov.theta <- cov(theta)
    
    #### posterior sd from covariance diagonal ####
    post.sd.vec <- sqrt(diag(cov.theta))
    post.sd <- matrix(post.sd.vec, d1, d2)
    
    #### credible intervals ####
    ci.lower <- apply(arr, c(1,2), quantile, probs = 0.025)
    ci.upper <- apply(arr, c(1,2), quantile, probs = 0.975)
    
    #### coverage ####
    coverage <- (true_value >= ci.lower) & (true_value <= ci.upper)
    
    #### RMSE ####
    rmse <- sqrt(Reduce("+", lapply(chain, function(x) (x - true_value)^2)) / n_iter)
    
  } else {
    
    vals <- unlist(chain)
    
    post.mean <- mean(vals)
    post.median <- median(vals)
    post.sd <- sd(vals)
    
    ci.lower <- quantile(vals, 0.025)
    ci.upper <- quantile(vals, 0.975)
    
    coverage <- (true_value >= ci.lower) & (true_value <= ci.upper)
    
    rmse <- sqrt(mean((vals - true_value)^2))
    
    cov.theta <- var(vals)
  }
  
  list(
    post.mean = post.mean,
    post.median = post.median,
    post.sd = post.sd,
    ci.lower = ci.lower,
    ci.upper = ci.upper,
    coverage = coverage,
    rmse = rmse,
    post.cov = cov.theta
  )
}



score_function <- function(W.obs.ord, Y.obs.ord, family) {
  q <- ncol(W.obs.ord)
  n <- nrow(W.obs.ord)
  scores <- matrix(NA, nrow = n, ncol = q)
  
  for (j in 1:q) {
    W_j <- W.obs.ord[, j]
    Y_j <- Y.obs.ord[, j]
    
    if (family[j] == "Gaussian") {
      scores[, j] <- Y_j - W_j
    } else if (family[j] == "Poisson") {
      scores[, j] <- Y_j - exp(W_j)
    } else if (family[j] == "Binomial") {
      p <- exp(W_j) / (1 + exp(W_j))
      scores[, j] <- Y_j - p
    } else if (family[j] == "Gamma") {
      alpha <- 1  # Can be customized or passed
      scores[, j] <- -Y_j * exp(W_j) + alpha
    } else if (family[j] == "Negative-Binomial") {
      mu <- exp(W_j)
      theta <- 1  # Fixed dispersion
      scores[, j] <- (Y_j - mu) / (mu + mu^2 / theta)
    } else if (family[j] == "Inverse-Gaussian") {
      mu <- exp(W_j)
      lambda <- 1
      scores[, j] <- (Y_j - mu) / (mu^3 / lambda)
    } else if (family[j] == "Tweedie") {
      mu <- exp(W_j)
      p <- 1.5
      scores[, j] <- Y_j - mu^p  # Approximate; subject to link function
    } else if (family[j] == "Multinomial") {
      p <- exp(W_j) / (1 + exp(W_j))
      scores[, j] <- Y_j - p
    } else {
      stop(paste("Unsupported family:", family[j]))
    }
  }
  
  return(vec(scores))
}


#### auto-tuning function for MH algorithm in warm up chain

tuning.update <- function(acc, att, tuning, nattempts = 50, lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.2) & these.update
  these.high   <- (acc.rate > 0.3) & these.update
  
  tuning[these.low]  <- tuning[these.low] * lower
  tuning[these.high] <- tuning[these.high] * higher
  
  acc[these.update] <- 0
  att[these.update] <- 0
  
  results <- list(acc = acc, att = att, tuning = tuning)
  return(results)
}
