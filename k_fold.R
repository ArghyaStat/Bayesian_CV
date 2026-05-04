library(doParallel)
library(foreach)
library(mcmcse)

set.seed(123)

# -----------------------------
# USER SETTINGS
# -----------------------------
n <- 2e2
p <- 1e1
K <- 10   # change to n for LOOCV
S <- 1e4
reps <- 100

sigma2 <- 1
tau2 <- 1e2

# -----------------------------
# DATA GENERATION (LOGISTIC)
# -----------------------------
sim_data <- function(n, p){
  X <- matrix(rnorm(n*p), n, p)
  beta <- rnorm(p)
  eta <- X %*% beta
  prob <- 1/(1+exp(-eta))
  y <- rbinom(n,1,prob)
  list(X=X,y=y)
}

# -----------------------------
# KAPPA FUNCTION
# -----------------------------
make_folds <- function(n,K){
  sample(rep(1:K, length.out=n))
}

# -----------------------------
# LOG-LIKELIHOOD (LOGISTIC)
# -----------------------------
log_like <- function(theta, X, y){
  eta <- X %*% theta
  y*eta - log(1+exp(eta))
}

# -----------------------------
# PRIOR
# -----------------------------
log_prior <- function(theta){
  sum(dnorm(theta,0,sqrt(tau2),log=TRUE))
}

log_post <- function(theta){
  sum(log_like(theta,X,y)) + log_prior(theta)
}

U <- function(theta) -log_post(theta)

grad <- function(theta){
  eta <- X %*% theta
  p_hat <- 1/(1+exp(-eta))
  t(X) %*% (p_hat - y) + theta/tau2
}

# Bronze target
U_b <- function(theta){
  ll <- log_like(theta,X,y)
  -( sum(ll)*(K-1)/K + log_prior(theta) )
}

# Zanella mixture target
U_z <- function(theta){
  ll <- log_like(theta,X,y)
  fold_ll <- tapply(ll, kappa, sum)
  -( sum(ll) + log(sum(exp(-fold_ll))) + log_prior(theta) )
}

# HMC sampler (copied from our paper on HMC for dummies)

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
# MAIN REPLICATION LOOP
# -----------------------------
cl <- makeCluster(10)
registerDoParallel(cl)

results <- foreach(r=1:reps, .combine=rbind,
                   .packages=c("mcmcse")) %dopar% {
                     
                     dat <- sim_data(n,p)
                     X <- dat$X
                     y <- dat$y
                     kappa <- make_folds(n,K)
                     
                     # -------------------------
                     # Gelfand samples
                     # -------------------------
                     log_post <- function(theta){
                       sum(log_like(theta,X,y)) + log_prior(theta)
                     }
                     
                     U <- function(theta) -log_post(theta)
                     
                     grad <- function(theta){
                       eta <- X %*% theta
                       p <- 1/(1+exp(-eta))
                       t(X) %*% (p - y) + theta/tau2
                     }
                     
                     hmc_g <- leapfrog_hmc(0.01,10,S,U,grad,rep(0,p),diag(p))
                     theta_g <- hmc_g$samples
                     acc_g <- hmc_g$acc
                     
                     hmc_b <- leapfrog_hmc(0.01,10,S,U_b,grad,rep(0,p),diag(p))
                     theta_b <- hmc_b$samples
                     acc_b <- hmc_b$acc
                     
                     hmc_z <- leapfrog_hmc(0.01,10,S,U_z,grad,rep(0,p),diag(p))
                     theta_z <- hmc_z$samples
                     acc_z <- hmc_z$acc
                     
                     # -------------------------
                     # Likelihood matrices
                     # -------------------------
                     lik_g <- exp(sapply(1:S, function(s) log_like(theta_g[s,],X,y)))
                     lik_b <- exp(sapply(1:S, function(s) log_like(theta_b[s,],X,y)))
                     lik_z <- exp(sapply(1:S, function(s) log_like(theta_z[s,],X,y)))
                     
                     lik_g <- t(lik_g); lik_b <- t(lik_b); lik_z <- t(lik_z)
                     
                     # -------------------------
                     # FUNCTION: compute ELPPD + variance
                     # -------------------------
                     compute_elppd <- function(lik, type){
                       
                       mu <- numeric(n)
                       num <- den <- matrix(NA,S,n)
                       
                       for(i in 1:n){
                         
                         fold <- which(kappa==kappa[i])
                         
                         if(type=="G"){
                           w <- 1 / apply(lik[,fold],1,prod)
                         }
                         
                         if(type=="B"){
                           w <- 1 / apply(lik[,fold],1,prod)
                         }
                         
                         if(type=="Z"){
                           denom <- rowSums(sapply(1:K, function(k){
                             idx <- which(kappa==k)
                             1 / apply(lik[,idx],1,prod)
                           }))
                           w <- (1/apply(lik[,fold],1,prod)) / denom
                         }
                         
                         num[,i] <- lik[,i]*w
                         den[,i] <- w
                         mu[i] <- sum(num[,i])/sum(den[,i])
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

colnames(results) <- c(
  "elpd_g","var_g","varS_g",
  "elpd_b","var_b","varS_b",
  "elpd_z","var_z","varS_z",
  "acc_g","acc_b","acc_z"
)

cat("\nAcceptance Rates (mean ± sd):\n")

cat(sprintf("Gelfand : %.3f ± %.3f\n",
            mean(results[,"acc_g"]),
            sd(results[,"acc_g"])))

cat(sprintf("Bronze  : %.3f ± %.3f\n",
            mean(results[,"acc_b"]),
            sd(results[,"acc_b"])))

cat(sprintf("Zanella : %.3f ± %.3f\n",
            mean(results[,"acc_z"]),
            sd(results[,"acc_z"])))
# -----------------------------
# PLOTS
# -----------------------------


# -----------------------------
# PLOTS (publication-ready)
# -----------------------------

cex_lab  <- 1.8
cex_axis <- 1.5
cex_main <- 1.6

# ---- PDF output (wide layout for 3 panels) ----
pdf("kfold_box_plots.pdf", width = 12, height = 4)

# ---- Layout and tight margins ----
par(mfrow = c(1,3),
    mar  = c(3.2, 3.8, 1.8, 0.8),   # bottom, left, top, right
    oma  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2, 0.7, 0),            # axis title, labels, line
    tcl  = -0.3,                    # tick length (short)
    cex.lab  = cex_lab,
    cex.axis = cex_axis,
    cex.main = cex_main)

# ---- Gelfand ----
boxplot(results[,"var_g"],
        main = "G (Gelfand)",
        ylab = "Estimated Variance",
        col = "gray90",
        border = "black")
abline(h = mean(results[,"var_g"]), col = "red", lwd = 2)

# ---- Bronze ----
boxplot(results[,"var_b"],
        main = "B (Bronze)",
        ylab = "",
        col = "gray90",
        border = "black")
abline(h = mean(results[,"var_b"]), col = "red", lwd = 2)

# ---- Zanella ----
boxplot(results[,"var_z"],
        main = "Z (Zanella)",
        ylab = "",
        col = "gray90",
        border = "black")
abline(h = mean(results[,"var_z"]), col = "red", lwd = 2)

dev.off()