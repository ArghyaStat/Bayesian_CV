update.W <- function(W.obs.ord, beta, Sigma, B, r, NNarray.obs,
                     Y.obs.ord, X.obs.ord, N.obs, log.like, family,
                     family.par = NULL,
                     link = NULL){
  
  q <- ncol(W.obs.ord)
  
  W.prior.mean.full <- X.obs.ord %*% beta
  
  post.W <- log.likelihood(W.obs.ord, Y.obs.ord, family,
                           family.par = family.par,
                           link = link)$like.total
  
  ## ---- draw all Gaussian directions once ----
  Z <- matrix(rnorm(N.obs*q), N.obs, q)
  
  ## ---- Vecchia draws ----
  nu <- vecchia_prior_draw(Z, NNarray.obs, B, r)
  
  for (j in seq_len(q)) {
    
    ind.minus <- setdiff(seq_len(q), j)
    
    if (length(ind.minus) > 0) {
      
      Sigma.jm <- Sigma[j, ind.minus, drop = FALSE]
      Sigma.mm <- Sigma[ind.minus, ind.minus, drop = FALSE]
      
      Sigma.mm.inv <- solve(Sigma.mm)
      
      cond.vec <- Sigma.mm.inv %*% Sigma[ind.minus, j, drop = FALSE]
      
      mu.j <- W.prior.mean.full[, j] +
        (W.obs.ord[, ind.minus, drop = FALSE] -
           W.prior.mean.full[, ind.minus, drop = FALSE]) %*%
        cond.vec
      
      Sigma.cond <- as.numeric(
        Sigma[j, j] -
          Sigma.jm %*% Sigma.mm.inv %*% t(Sigma.jm)
      )
      
    } else {
      
      mu.j <- W.prior.mean.full[, j]
      Sigma.cond <- Sigma[j, j]
      
    }
    
    ## ESS direction
    W.prior.j <- sqrt(Sigma.cond) * nu[,j]
    
    gamma <- runif(1, 0, 2*pi)
    gamma_min <- gamma - 2*pi
    gamma_max <- gamma
    
    logy <- post.W + log(runif(1))
    
    while(TRUE){
      
      can.W.obs.ord <- W.obs.ord
      
      can.W.obs.ord[, j] <-
        mu.j +
        (W.obs.ord[, j] - mu.j) * cos(gamma) +
        W.prior.j * sin(gamma)
      
      post.can.W <- log.likelihood(can.W.obs.ord, Y.obs.ord, family,
                                   family.par = family.par,
                                   link = link)$like.total
      
      log.r.W <- post.can.W - post.W
      
      if (log(runif(1)) < log.r.W) {
        
        W.obs.ord <- can.W.obs.ord
        post.W <- post.can.W
        break
        
      } else {
        
        if (gamma < 0) gamma_min <- gamma else gamma_max <- gamma
        gamma <- runif(1, gamma_min, gamma_max)
        
      }
    }
  }
  
  W.obs.ord
}



update.phi <- function(phi, beta, chol.Sigma, nu,
                       W.obs.ord, distobs.nn, distall.nn,
                       NNarray.obs, NNarray.all,
                       coeff.all, acc.phi, tuning.phi, b_phi){
  
  N.obs <- nrow(W.obs.ord)
  
  # extract coeff.obs from coeff.all (current state)
  coeff.obs <- list(
    B = coeff.all$B[1:N.obs],
    r = coeff.all$r[1:N.obs]
  )
  
  # current posterior
  post.phi <- vecchia_ll(W.ord = W.obs.ord, 
                         M.ord = X.obs.ord %*% beta, 
                         NNarray = NNarray.obs, 
                         coeff = coeff.obs, 
                         chol.V = chol.Sigma, 
                         log = TRUE) + dunif(phi, 0, b_phi, log = TRUE)
  
  # proposal
  can.phi <- rtruncnorm(1, a = 0, b = b_phi, 
                        mean = phi, sd = sqrt(tuning.phi))
  
 
  can.coeff.all <- vecchia_coeff(distall.nn, NNarray.all, can.phi, nu)
  
  # extract obs part
  can.coeff.obs <- list(
    B = can.coeff.all$B[1:N.obs],
    r = can.coeff.all$r[1:N.obs]
  )
  
  # candidate posterior
  post.can.phi <- vecchia_ll(W.ord = W.obs.ord, 
                             M.ord = X.obs.ord %*% beta, 
                             NNarray = NNarray.obs, 
                             coeff = can.coeff.obs, 
                             chol.V = chol.Sigma, 
                             log = TRUE) + dunif(can.phi, 0, b_phi, log = TRUE)
  
  # proposal densities
  q_curr_given_can <- dtruncnorm(phi, a = 0, b = b_phi, 
                                 mean = can.phi, sd = sqrt(tuning.phi))
  q_can_given_curr <- dtruncnorm(can.phi, a = 0, b = b_phi, 
                                 mean = phi, sd = sqrt(tuning.phi))
  
  log.r.phi <- (post.can.phi - post.phi) + 
    log(q_curr_given_can) - log(q_can_given_curr)
  
  # accept/reject
  if(log(runif(1)) < log.r.phi){
    phi <- can.phi
    coeff.all <- can.coeff.all  
    acc.phi <- acc.phi + 1
  }
  
  return(list(
    phi = phi,
    coeff.all = coeff.all,
    acc.phi = acc.phi
  ))
}



update.Sigma <- function(UW,
                         M.prior, chol.prec.V.prior,
                         M.tilde.part, chol.V.tilde,
                         S.prior, df.prior, N.obs,
                         family,
                         fixed_variances) {
  
  q <- ncol(UW)
  
  rss.W.obs  <- crossprod(UW)
  
  rss.beta.prior <- crossprod(chol.prec.V.prior %*% M.prior)
  rss.beta.post  <- crossprod(chol.V.tilde %*% M.tilde.part)
  
  S.tilde <- S.prior + rss.W.obs + rss.beta.prior - rss.beta.post
  df.tilde <- df.prior + N.obs
  
  Sigma.star <- riwish(df.tilde, S.tilde)
  
  d.star <- sqrt(diag(Sigma.star))
  R.star <- diag(1/d.star) %*% Sigma.star %*% diag(1/d.star)
  
  D.new <- diag(q)
  
  for (j in 1:q) {
    if (family[j] == "Binomial") {
      D.new[j,j] <- sqrt(fixed_variances[j])
    } else {
      D.new[j,j] <- d.star[j]
    }
  }
  
  Sigma <- D.new %*% R.star %*% D.new
  
  chol.Sigma <- chol(Sigma)
  
  return(list(
    Sigma = Sigma,
    chol.Sigma = chol.Sigma,
    chol.Sigma.inv = chol(chol2inv(chol.Sigma))
  ))
}

# beta update

update.beta <- function(M.tilde.part, chol.Sigma, chol.V.tilde){

  beta.samp <- (matrix(rnorm(p*q, 0, 1), nrow = p , ncol = q) %*% chol.Sigma)

  beta.temp <- (chol.V.tilde %*% M.tilde.part + beta.samp)

  beta <- crossprod(chol.V.tilde, beta.temp)
  
  beta
  
}