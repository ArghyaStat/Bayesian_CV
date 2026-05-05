rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

source("loocv_functions.R")

set.seed(1)

library(mcmcse)

#  Data Simulation (Model: simple linear regression)

# True parameters theta_1 = 1, theta_2 = 2
# Data size = 50

n <- 50
x <- rnorm(n)
y <- 1 + 2*x + rnorm(n)
sigma2 <- 1
tau2 <- 100

# Gelfand sampler via HMC 

# MCMC sample size S = 5e3 (user's choice)
S <- 5e3

hmc_gelfand <- leapfrog_hmc(
  epsilon = 0.25,
  L = 20,
  niters = S,
  U = U_post,
  grad = grad_U_post,
  x_init = c(1, 2),
  M = diag(2)
)

acc_gelfand <- hmc_gelfand$acc # Optimal acc = 0.65
cat(sprintf("HMC Gelfand acceptance rate: %.4f\n", acc_gelfand))
theta_g <- hmc_gelfand$samples


# Zanella (HMC on q_mix)
hmc_zanella <- leapfrog_hmc(
  epsilon = 0.24, L = 20, niters = S,
  U = U_mix, grad = grad_U_mix,
  x_init = c(1,2), M = diag(2)
)

acc_zanella <- hmc_zanella$acc
cat(sprintf("HMC Zanella acceptance rate: %.4f\n", acc_zanella))
theta_z <- hmc_zanella$samples


# Log-likelihood matrices (FIXED) of dimension S \times n for Gelfand and Zanella
# S = no.of MCMC samples, n = data size

loglik_g <- loglik_z <- matrix(NA, S, n)

for(s in 1:S){
  loglik_g[s,] <- log_like(theta_g[s,])
  loglik_z[s,] <- log_like(theta_z[s,])
}

lik_g <- exp(loglik_g)
lik_z <- exp(loglik_z)

# LOOCV estimators

# mu_hat and ELPD from Gelfand ----

mu_g <- sapply(1:n, function(i) 1/mean(1/lik_g[,i]))
elpd_g <- sum(log(mu_g))


## Important sampling SE

# Gelfand
num_g <- matrix(1, S, n)
den_g <- 1/lik_g
grad_g <- c(rep(1, n), -mu_g)

vec_g <- cbind(num_g, den_g)

foo_g <- mcse.multi(vec_g, r = 1, adjust = FALSE, size = 100)
Omega_hat_g <- foo_g$cov
se_g <- as.numeric(t(grad_g) %*% Omega_hat_g %*% grad_g / S)

# mu_hat and ELPD from Zanella (self-normalized mixture) ----

mu_z <- numeric(n)
num_z <- den_z <- matrix(NA, S, n)

for(i in 1:n){
  w <- exp(-loglik_z[,i])
  w <- w / sum(w)
  num_z[,i] <- lik_z[,i] * w
  den_z[,i] <- w
  mu_z[i] <- sum(num_z[,i])/sum(w)
}

elpd_z <- sum(log(mu_z))


# Zanella (self-normalized)
g1 <- 1/colMeans(num_z)
g2 <- -1/colMeans(den_z)
grad_z <- c(g1, g2)

foo_z <- mcse.multi(cbind(num_z, den_z), r = 1, adjust = FALSE, size = 100)
Omega_hat_z <- foo_z$cov
se_z <- as.numeric(t(grad_z) %*% Omega_hat_z %*% grad_z / S)

# PSIS

mu_p <- numeric(n)
num_p <- den_p <- matrix(NA, S, n)

for(i in 1:n){
  w <- psis_smooth(log_ratios = -loglik_g[,i], tail_frac = 0.2)
  num_p[,i] <- lik_g[,i] * w
  den_p[,i] <- w
  mu_p[i] <- sum(num_p[,i])/sum(den_p[,i])
}

elpd_p <- sum(log(mu_p))

grad_p <- c(1/colMeans(num_p), -1/colMeans(den_p))
foo_p <- mcse.multi(cbind(num_p, den_p), adjust = FALSE, r = 1)
se_p <- as.numeric(t(grad_p) %*% foo_p$cov %*% grad_p / S)



reps <- 1e2

elpd_rep_g <- numeric(reps)
elpd_rep_z <- numeric(reps)
elpd_rep_p <- numeric(reps)

start_time <- Sys.time()  # Start time using Sys.time()
for (r in 1:reps) {
  
  if(r %% (reps / 10) == 0) {
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("Main Phase Progress: %.0f%%, Elapsed time: %.2f seconds\n", 
                100 * (r / reps), elapsed_time))
  }

  # Gelfand
  
  theta_g <- leapfrog_hmc(
    epsilon = 0.25,
    L = 20,
    niters = S,
    U = U_post,
    grad = grad_U_post,
    x_init = c(1, 2),
    M = diag(2)
  )$samples
  
  lik_g <- matrix(NA, S, n)
  for (s in 1:S) {
    lik_g[s, ] <- exp(log_like(theta_g[s, ]))
  }
  
  loocv_g <- numeric(n)
  for (i in 1:n) {
    loocv_g[i] <- 1 / mean(1 / lik_g[, i])
  }
  
  elpd_rep_g[r] <- sum(log(loocv_g))
  
  
  # Zanella (q_mix via HMC)
  
  theta_z <- leapfrog_hmc(
    epsilon = 0.24, L = 20, niters = S,
    U = U_mix, grad = grad_U_mix,
    x_init = c(1,2), M = diag(2)
  )$samples
  
  lik_z <- matrix(NA, S, n)
  for (s in 1:S) {
    lik_z[s, ] <- exp(log_like(theta_z[s, ]))
  }
  
  loocv_z <- numeric(n)
  
  for (i in 1:n) {
    w <- 1 / lik_z[, i]
    w <- w / sum(w)   # self-normalized mixture weights
    loocv_z[i] <- sum(lik_z[, i] * w)/sum(w)
  }
  
  elpd_rep_z[r] <- sum(log(loocv_z))
  
  loocv_p <- sapply(1:n, function(i){
    w <- psis_smooth(log_ratios = -log(lik_g[,i]), tail_frac = 0.2)
    sum(lik_g[,i]*w)/sum(w)
  })
  
  elpd_rep_p[r] <- sum(log(loocv_p))

}
end_time <- Sys.time()
total_elapsed_time <- difftime(end_time, start_time, units = "secs")[[1]] 

elpd_outputs <- list("elpd_gelfand_is" = elpd_g,
                     "se_elpd_gelfand_is" = se_g,
                     "elpd_zanella_is" = elpd_z,
                     "se_elpd_zanella_is" = se_z,
                     "elpd_psis_is" = elpd_p,
                     "se_elpd_psis_is" = se_p,
                     "elpd_gelfand_rep" = elpd_rep_g,
                     "elpd_zanella_rep" = elpd_rep_z,
                     "elpd_psis_rep" = elpd_rep_p)

save(elpd_outputs, file = "elpd_outputs.RData")

# ---- Save Posterior Samples ----
# save(theta_g, acc_gelfand, file = "HMC_gelfand.RData")
# load("HMC_gelfand.RData") 


# ---- Save Posterior Samples ----
# save(theta_z, acc_zanella, file = "HMC_Zanella.RData")
# load("HMC_Zanella.RData") 


acc_g <- theta_g$acc
acc_b <- theta_b$acc
acc_z <- theta_z$acc
