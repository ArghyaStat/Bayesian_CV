
rm(list = ls())

library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)

# -----------------------------
# SETTINGS
# -----------------------------
n <- 50
p <- 10
K <- 5
S <- 1e4
sigma2 <- 1e2

# -----------------------------
# DATA (FIXED)
# -----------------------------
sim_data <- function(n, p){
  X <- matrix(rnorm(n*p), n, p)
  beta <- rnorm(p)
  eta <- X %*% beta
  prob <- 1/(1+exp(-eta))
  y <- rbinom(n,1,prob)
  list(X=X,y=y)
}

dat <- sim_data(n,p)
X <- dat$X
y <- dat$y

kappa <- sample(rep(1:K, length.out=n))

# -----------------------------
# HELPERS
# -----------------------------
log_sum_exp <- function(x){
  m <- max(x)
  m + log(sum(exp(x - m)))
}

log_like <- function(theta, X, y){
  eta <- X %*% theta
  y*eta - log(1+exp(eta))
}

log_prior <- function(theta){
  sum(dnorm(theta,0,sqrt(sigma2),log=TRUE))
}

# -----------------------------
# GRADIENTS
# -----------------------------
grad_post <- function(theta){
  eta <- X %*% theta
  p_hat <- 1/(1+exp(-eta))
  -( t(X) %*% (y - p_hat) + theta/sigma2 )
}

grad_b <- function(theta){
  eta <- X %*% theta
  p_hat <- 1/(1+exp(-eta))
  -( (K-1)/K * t(X) %*% (y - p_hat) + theta/sigma2 )
}

grad_z <- function(theta){
  
  eta <- X %*% theta
  p_hat <- 1/(1+exp(-eta))
  r <- y - p_hat
  
  grad_ll <- t(X) %*% r
  
  ll_vec <- log_like(theta,X,y)
  fold_ll <- tapply(ll_vec, kappa, sum)
  
  log_w <- -fold_ll
  log_w <- log_w - log_sum_exp(log_w)
  w <- exp(log_w)
  
  grad_mix <- 0
  for(k in 1:K){
    idx <- which(kappa==k)
    grad_k <- t(X[idx,,drop=FALSE]) %*% r[idx]
    grad_mix <- grad_mix - w[k] * grad_k
  }
  
  -(grad_ll + grad_mix + theta/sigma2)
}

# -----------------------------
# TARGETS
# -----------------------------
U_post <- function(theta){
  -(sum(log_like(theta,X,y)) + log_prior(theta))
}

U_b <- function(theta){
  ll <- sum(log_like(theta,X,y))
  -( (K-1)/K * ll + log_prior(theta) )
}

U_z <- function(theta){
  ll_vec <- log_like(theta,X,y)
  fold_ll <- tapply(ll_vec, kappa, sum)
  -( sum(ll_vec) + log_sum_exp(-fold_ll) + log_prior(theta) )
}

# -----------------------------
# HMC
# -----------------------------
leapfrog_hmc <- function(epsilon, L, niters, U, grad, x_init, M) {
  d <- length(x_init)
  M_inv <- chol2inv(chol(M))
  x <- matrix(NA, niters, d)
  x[1,] <- x_init
  
  chol_M <- chol(M)
  p <- matrix(rnorm(niters*d), niters, d) %*% chol_M
  acc <- 0
  
  for(k in 2:niters){
    x_prop <- x[k-1,]
    p_prop <- p[k,]
    
    g <- grad(x_prop)
    p_prop <- p_prop - 0.5*epsilon*g
    
    for(l in 1:L){
      x_prop <- x_prop + epsilon*(M_inv %*% p_prop)
      g <- grad(x_prop)
      if(l != L) p_prop <- p_prop - epsilon*g
    }
    p_prop <- p_prop - 0.5*epsilon*g
    
    H_curr <- U(x[k-1,]) + 0.5*t(p[k,]) %*% (M_inv %*% p[k,])
    H_prop <- U(x_prop) + 0.5*t(p_prop) %*% (M_inv %*% p_prop)
    
    if(log(runif(1)) < -(H_prop - H_curr)){
      x[k,] <- x_prop
      acc <- acc + 1
    } else {
      x[k,] <- x[k-1,]
    }
  }
  list(samples = x, acc = acc/niters)
}


reps <- 100

cl <- makeCluster(10)
registerDoParallel(cl)

results <- foreach(r = 1:reps, .combine = rbind,
                   .packages = c("mcmcse")) %dopar% {
                     
                     # -----------------------------
                     # HMC SAMPLING (PER REP)
                     # -----------------------------
                     x_init <- glm(y ~ X - 1, family="binomial")$coef
                     
                     theta_g <- leapfrog_hmc(0.35,10,S,U_post,grad_post,x_init,diag(p))
                     theta_b <- leapfrog_hmc(0.35,10,S,U_b,grad_b,x_init,diag(p))
                     theta_z <- leapfrog_hmc(0.15,10,S,U_z,grad_z,x_init,diag(p))
                     
                     samples_g <- theta_g$samples
                     samples_b <- theta_b$samples
                     samples_z <- theta_z$samples
                     
                     acc_g <- theta_g$acc
                     acc_b <- theta_b$acc
                     acc_z <- theta_z$acc
                     
                     # -----------------------------
                     # LIKELIHOODS
                     # -----------------------------
                     lik_g <- t(exp(sapply(1:S, function(s) log_like(samples_g[s,],X,y))))
                     lik_b <- t(exp(sapply(1:S, function(s) log_like(samples_b[s,],X,y))))
                     lik_z <- t(exp(sapply(1:S, function(s) log_like(samples_z[s,],X,y))))
                     
                     # -----------------------------
                     # FOLD LIKELIHOODS
                     # -----------------------------
                     compute_fold_lik <- function(lik){
                       sapply(1:K, function(k){
                         idx <- which(kappa==k)
                         apply(lik[,idx,drop=FALSE],1,prod)
                       })
                     }
                     
                     fold_lik_z <- compute_fold_lik(lik_z)
                     
                     # -----------------------------
                     # Gelfand
                     # -----------------------------
                     mu_g <- numeric(n)
                     num_g <- matrix(1, S, n)
                     den_g <- 1 / lik_g
                     
                     for(i in 1:n){
                       mu_g[i] <- 1 / mean(den_g[,i])
                     }
                     
                     elpd_g <- sum(log(mu_g))
                     
                     grad_g <- c(rep(1,n), -mu_g)
                     Omega_g <- mcse.multi(cbind(num_g, den_g), r=1)$cov
                     var_g <- as.numeric(t(grad_g) %*% Omega_g %*% grad_g / S)
                     varS_g <- var(rowSums(log(num_g / den_g))) / S
                     
                     # -----------------------------
                     # Bronze
                     # -----------------------------
                     mu_b <- numeric(n)
                     num_b <- matrix(1, S, n)
                     den_b <- 1 / lik_b
                     
                     for(i in 1:n){
                       mu_b[i] <- 1 / mean(den_b[,i])
                     }
                     
                     elpd_b <- sum(log(mu_b))
                     
                     grad_b <- c(rep(1,n), -mu_b)
                     Omega_b <- mcse.multi(cbind(num_b, den_b), r=1)$cov
                     var_b <- as.numeric(t(grad_b) %*% Omega_b %*% grad_b / S)
                     varS_b <- var(rowSums(log(num_b / den_b))) / S
                     
                     # -----------------------------
                     # Zanella
                     # -----------------------------
                     mu_z <- numeric(n)
                     num_z <- den_z <- matrix(NA, S, n)
                     
                     denom_all <- rowSums(1 / fold_lik_z)
                     
                     for(i in 1:n){
                       
                       k_i <- kappa[i]
                       
                       w <- (1 / fold_lik_z[,k_i]) / denom_all
                       
                       num_z[,i] <- lik_z[,i] * w
                       den_z[,i] <- w
                       
                       mu_z[i] <- sum(num_z[,i]) / sum(den_z[,i])
                     }
                     
                     elpd_z <- sum(log(mu_z))
                     
                     grad_z_vec <- c(1/colMeans(num_z), -1/colMeans(den_z))
                     Omega_z <- mcse.multi(cbind(num_z, den_z), r=1)$cov
                     var_z <- as.numeric(t(grad_z_vec) %*% Omega_z %*% grad_z_vec / S)
                     varS_z <- var(rowSums(log(num_z / den_z))) / S
                     
                     # -----------------------------
                     # RETURN
                     # -----------------------------
                     c(
                       elpd_g, var_g, varS_g,
                       elpd_b, var_b, varS_b,
                       elpd_z, var_z, varS_z,
                       acc_g, acc_b, acc_z
                     )
                   }
# -----------------------------
# COLUMN NAMES
# -----------------------------
colnames(results) <- c(
  "elpd_g","var_g","varS_g",
  "elpd_b","var_b","varS_b",
  "elpd_z","var_z","varS_z",
  "acc_g","acc_b","acc_z"
)

# -----------------------------
# SAVE
# -----------------------------
save(results, file = paste0(K, "_fold_elppd_logistic.RData"))