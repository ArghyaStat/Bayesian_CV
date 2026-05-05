rm(list = ls())

mydir <- this.path::here()
setwd(mydir)


library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)

# -----------------------------
# USER SETTINGS
# -----------------------------
n <- 50
p <- 10
K <- 10
S <- 1e4
reps <- 100

sigma2 <- 1e2

# -----------------------------
# DATA
# -----------------------------
sim_data <- function(n, p){
  X <- matrix(rnorm(n*p), n, p)
  beta <- rnorm(p)
  eta <- X %*% beta
  prob <- 1/(1+exp(-eta))
  y <- rbinom(n,1,prob)
  list(X=X,y=y)
}

make_folds <- function(n,K){
  sample(rep(1:K, length.out=n))
}

# -----------------------------
# MODEL
# -----------------------------
log_like <- function(theta, X, y){
  eta <- X %*% theta
  y*eta - log(1+exp(eta))
}

log_prior <- function(theta){
  sum(dnorm(theta,0,sqrt(sigma2),log=TRUE))
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

# -----------------------------
# PARALLEL
# -----------------------------
cl <- makeCluster(K)
registerDoParallel(cl)

results <- foreach(r=1:reps, .combine=rbind,
                   .packages=c("mcmcse")) %dopar% {
                     
                     dat <- sim_data(n,p)
                     X <- dat$X
                     y <- dat$y
                     kappa <- make_folds(n,K)
                     
                     # -------------------------
                     # Gradients
                     # -------------------------
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
                       
                       w <- exp(-fold_ll)
                       w <- w / sum(w)
                       
                       grad_mix <- 0
                       for(k in 1:K){
                         idx <- which(kappa==k)
                         grad_k <- t(X[idx,,drop=FALSE]) %*% r[idx]
                         grad_mix <- grad_mix - w[k]*grad_k
                       }
                       
                       -(grad_ll + grad_mix + theta/sigma2)
                     }
                     
                     # -------------------------
                     # Targets
                     # -------------------------
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
                       -( sum(ll_vec) + log(sum(exp(-fold_ll))) + log_prior(theta) )
                     }
                     
                     # -------------------------
                     # HMC
                     # -------------------------
                     hmc_g <- leapfrog_hmc(epsilon = 1e-2, 
                                           L = 10, 
                                           niters = S, 
                                           U = U_post, 
                                           grad = grad_post, 
                                           x_init = rep(0,p), 
                                           M = diag(p))
                       
                       
                     theta_g <- hmc_g$samples
                     acc_g <- hmc_g$acc
                     
                     hmc_b <- leapfrog_hmc(epsilon = 1e-2, 
                                           L = 10, 
                                           niters = S, 
                                           U = U_b, 
                                           grad = grad_b, 
                                           x_init = rep(0,p), 
                                           M = diag(p))
                     theta_b <- hmc_b$samples
                     acc_b <- hmc_b$acc
                     
                     hmc_z <- leapfrog_hmc(epsilon = 1e-2, 
                                           L = 10, 
                                           niters = S, 
                                           U = U_z, 
                                           grad = grad_z, 
                                           x_init = rep(0,p), 
                                           M = diag(p))
                     theta_z <- hmc_z$samples
                     acc_z <- hmc_z$acc
                     
                     # -------------------------
                     # Likelihoods
                     # -------------------------
                     lik_g <- t(exp(sapply(1:S, function(s) log_like(theta_g[s,],X,y))))
                     lik_b <- t(exp(sapply(1:S, function(s) log_like(theta_b[s,],X,y))))
                     lik_z <- t(exp(sapply(1:S, function(s) log_like(theta_z[s,],X,y))))
                     
                     # -------------------------
                     # ELPPD
                     # -------------------------
                     compute_elppd <- function(lik, type){
                       
                       mu <- numeric(n)
                       num <- den <- matrix(NA,S,n)
                       
                       for(i in 1:n){
                         fold <- which(kappa==kappa[i])
                         lik_fold <- apply(lik[,fold,drop=FALSE],1,prod)
                         
                         if(type %in% c("G","B")){
                           inv_lik <- 1/lik_fold
                           mu[i] <- 1/mean(inv_lik)
                           
                           num[,i] <- 1
                           den[,i] <- inv_lik
                         }
                         
                         if(type=="Z"){
                           denom_all <- rowSums(sapply(1:K,function(k){
                             idx <- which(kappa==k)
                             1/apply(lik[,idx,drop=FALSE],1,prod)
                           }))
                           
                           w <- (1/lik_fold)/denom_all
                           w <- w/sum(w)
                           
                           num[,i] <- lik[,i]*w
                           den[,i] <- w
                           
                           mu[i] <- sum(num[,i])/sum(den[,i])
                         }
                       }
                       
                       elppd <- sum(log(mu))
                       
                       grad <- c(1/colMeans(num), -1/colMeans(den))
                       foo <- mcse.multi(cbind(num,den), r=1, adjust=FALSE)
                       var_clt <- as.numeric(t(grad) %*% foo$cov %*% grad / S)
                       
                       var_sample <- var(rowSums(log(num/den)))/S
                       
                       c(elppd, var_clt, var_sample)
                     }
                     
                     res_g <- compute_elppd(lik_g,"G")
                     res_b <- compute_elppd(lik_b,"B")
                     res_z <- compute_elppd(lik_z,"Z")
                     
                     c(res_g, res_b, res_z, acc_g, acc_b, acc_z)
                   }

stopCluster(cl)

# -----------------------------
# OUTPUT
# -----------------------------
colnames(results) <- c(
  "elpd_g","var_g","varS_g",
  "elpd_b","var_b","varS_b",
  "elpd_z","var_z","varS_z",
  "acc_g","acc_b","acc_z"
)

save(results, file="k_fold_elppd_logistic.RData")