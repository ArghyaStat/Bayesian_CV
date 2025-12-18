rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

set.seed(123)

library(mcmcse)


# Simulation from data (linear regression)

n <- 50
x <- rnorm(n)
theta.true <- c(1, 2)
sigma <- 1
y <- theta.true[1] + theta.true[2] * x + rnorm(n, sd = sigma)

# Define log likelihood and log posterior

log_post <- function(theta){
  beta0 <- theta[1]
  beta1 <- theta[2]
  ll <- sum(dnorm(y, beta0 + beta1 * x, sd = sigma, log = TRUE))
  lp <- sum(dnorm(theta, 0, sqrt(100), log = TRUE))
  ll + lp
}

run_mcmc <- function(S, theta.init, proposal.sd){
  
  chain <- matrix(NA, S, 2)
  theta <- theta.init
  lp <- log_post(theta)
  
  acc <- 0
  
  for(s in 1:S){
    theta.prop <- theta + rnorm(2, 0, proposal.sd)
    lp.prop <- log_post(theta.prop)
    
    if(log(runif(1)) < lp.prop - lp){
      theta <- theta.prop
      lp <- lp.prop
      acc <- acc + 1
    }
    chain[s,] <- theta
  }
  cat(sprintf("MH acceptance rate: %.4f\n", (acc/S)))
  return(chain)
}

### Running MCMC

S <- 5e4
theta <- run_mcmc(S, theta.init = theta.true, proposal.sd = c(0.2, 0.2))

### Save likelihoods for S MCMC samples for n data points

lik <- matrix(NA, S, n)
for(i in 1:n){
  lik[,i] <- dnorm(y[i], theta[,1] + theta[,2] * x[i], sd = sigma)
}


# LOOCV estimators (return loo, num, den)

# ---- Gelfand ----
loo_gelfand <- function(lik){
  list(
    loo = apply(lik, 2, function(z) 1 / mean(1 / z)),
    num = matrix(1, nrow(lik), ncol(lik)),
    den = 1 / lik
  )
}

# ---- Zanella ----

loo_zanella <- function(lik){
  
  inv <- 1 / lik
  w <- inv / rowSums(inv)
  list(
    loo = colSums(lik * w)/colSums(w),
    num = lik * w,
    den = w
  )
}


# Hard-coded PSIS (from pseudocode) with help of ChatGPT

psis_smooth_1d <- function(r, tail_frac = 0.3){
  
  S <- length(r)
  M <- ceiling(tail_frac * S)
  
  ord <- order(r)
  r_sorted <- r[ord]
  
  tail <- r_sorted[(S - M + 1):S]
  u <- min(tail)
  y <- tail - u
  
  # Moment-based GPD fit
  k_hat <- 0.5 * (mean(y)^2 / var(y) - 1)
  sigma_hat <- mean(y) * (1 + k_hat)
  
  p <- (1:M - 0.5) / M
  tail_smooth <- u + (sigma_hat / k_hat) * ((1 - p)^(-k_hat) - 1)
  
  r_sorted[(S - M + 1):S] <- tail_smooth
  r_smooth <- r_sorted[order(ord)]
  
  # Truncation
  r_smooth <- pmin(r_smooth, sqrt(S) * mean(r_smooth))
  r_smooth
}

loo_psis <- function(lik){
  S <- nrow(lik)
  n <- ncol(lik)
  
  num <- matrix(NA, S, n)
  den <- matrix(NA, S, n)
  loo <- numeric(n)
  
  for(i in 1:n){
    r_raw <- 1 / lik[,i]
    r_smooth <- psis_smooth_1d(r_raw)
    w <- r_smooth / sum(r_smooth)
    
    num[,i] <- lik[,i] * w
    den[,i] <- w
    loo[i] <- sum(num[,i])/sum(den[,i])
  }
  list(loo = loo, num = num, den = den)
}


# Importance sampling variance via delta method with MCSE for Omega_hat

is_variance <- function(estimator){
  
  num <- estimator$num
  den <- estimator$den
  
  u_bar_inv <- S / colSums(num)
  w_bar_inv <- -S / colSums(den)
  grad <- c(u_bar_inv, w_bar_inv)
  
  z <- cbind(num, den)
  Omega.hat <- mcse.multi(z, r = 1, adjust = FALSE, size = 100)$cov
  
  as.numeric(t(grad) %*%  Omega.hat %*% grad / S)
}

##  Comparison of variance with replicated data

### Here we obtain ELPD from different estimators 
# We draw different sets of independent replications and compute variances

Rmax <- 5e2
replications <- c(1e2, 2e2, 5e2)

elpd_g <- elpd_z <- elpd_p <- numeric(Rmax)

for(r in 1:Rmax){
  
  k <- sample(S, replace = TRUE)
  
  elpd_g[r] <- sum(log(loo_gelfand(lik[k, ])$loo))
  elpd_z[r] <- sum(log(loo_zanella(lik[k, ])$loo))
  elpd_p[r] <- sum(log(loo_psis(lik[k, ])$loo))
  
}

mc_var <- sapply(replications, function(R){
  c(
    var(elpd_g[1:R]),
    var(elpd_z[1:R]),
    var(elpd_p[1:R])
  )
})


# Importace sampling variances for different methods


var_g <- is_variance(loo_gelfand(lik)) # Gelfand
var_z <- is_variance(loo_zanella(lik)) # Zanella
var_p <- is_variance(loo_psis(lik)) # PSIS


# plot


pdf("mc_variance_loocv.pdf", width = 6, height = 6)

par(
  mar = c(4, 5, 1, 1), oma = c(0,0,0,0)
)

plot(
  replications, mc_var[1,],
  type = "b", pch = 19, col = "black",
  ylim = range(c(mc_var, var_g, var_z, var_p)),
  xlab = "Number of replications",
  ylab = "Monte Carlo variance",
)

lines(replications, mc_var[2,], type = "b", pch = 17, col = "blue")
lines(replications, mc_var[3,], type = "b", pch = 15, col = "red")

abline(h = var_g, lty = 3, col = "black", lwd = 2)
abline(h = var_z, lty = 3, col = "blue",  lwd = 2)
abline(h = var_p, lty = 3, col = "red",   lwd = 2)

legend(
  "topright",
  legend = c(
    "MC Gelfand", "MC Zanella", "MC PSIS",
    "IS Gelfand", "IS Z", "IS PSIS"
  ),
  col = c("black","blue","red","black","blue","red"),
  pch = c(19,17,15,NA,NA,NA),
  lty = c(1,1,1,2,2,2),
  bty = "n"
)

dev.off()