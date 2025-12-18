set.seed(1)

############################################################
# 1. Simulate data
############################################################
n <- 50
x <- rnorm(n)
y <- 1 + 2 * x + rnorm(n, sd = 1)

sigma2 <- 1         # known variance of observation noise
tau2   <- 100       # known prior variance for beta0, beta1

############################################################
# 2. Log posterior for (beta0, beta1)
############################################################

log_posterior <- function(params, x, y, sigma2, tau2){
  beta0 <- params[1]
  beta1 <- params[2]
  
  # Priors: N(0, tau2)
  lp <- dnorm(beta0, 0, sqrt(tau2), log = TRUE) +
    dnorm(beta1, 0, sqrt(tau2), log = TRUE)
  
  # Likelihood
  mu <- beta0 + beta1 * x
  lp <- lp + sum(dnorm(y, mu, sqrt(sigma2), log = TRUE))
  
  return(lp)
}

############################################################
# 3. Bivariate random-walk Metropolis–Hastings
############################################################

run_mcmc <- function(x, y, sigma2, tau2, n_iter = 10000,
                     proposal_sd = c(0.1, 0.1)){
  
  chain <- matrix(NA, n_iter, 2)
  colnames(chain) <- c("beta0", "beta1")
  
  current <- c(1, 2)
  current_lp <- log_posterior(current, x, y, sigma2, tau2)
  
  for (t in 1:n_iter){
    proposal <- current + rnorm(2, 0, proposal_sd)
    proposal_lp <- log_posterior(proposal, x, y, sigma2, tau2)
    
    log_alpha <- proposal_lp - current_lp
    
    if (log(runif(1)) < log_alpha){
      current <- proposal
      current_lp <- proposal_lp
    }
    
    chain[t, ] <- current
  }
  
  chain
}

############################################################
# 4. Run MCMC
############################################################

mcmc <- run_mcmc(x, y, sigma2, tau2, n_iter = 1e5)

beta0 <- mcmc[,1]
beta1 <- mcmc[,2]
S <- nrow(mcmc)

############################################################
# 5. Gelfand LOOCV
############################################################

loocv <- numeric(n)
num <- matrix(1, nrow = S, ncol = n)
den <- matrix(0, nrow = S, ncol = n)
for (i in 1:n){
  mu_i <- beta0 + beta1 * x[i]
  
  # p(y_i | theta_s)
  lik_i <- dnorm(y[i], mean = mu_i, sd = sqrt(sigma2))
  
  # Gelfand estimator
  loocv[i] <- 1 / mean(1 / lik_i)
  
  den[,i] <- 1/lik_i
}

elpd <- sum(log(loocv))

############################################################
# 6. Output
############################################################

print(head(loocv))
cat("\nELPD (Gelfand LOOCV):", elpd, "\n")

############################################################
# 7. SE Calculation
############################################################
vec1 <- cbind(num,den)
grad.g1g2 <- c(rep(1, n), -loocv)

foo <- mcse.multi(vec1, r = 1, adjust = FALSE, size = 100)

t(grad.g1g2) %*% foo$cov %*% grad.g1g2/S

############################################################
# 7. With replication
############################################################

reps <- 1000
elpd_rep <- numeric(reps)
for(r in 1:reps)
{
  print(r)
  mcmc <- run_mcmc(x, y, sigma2, tau2, n_iter = 1e3)
  
  beta0 <- mcmc[,1]
  beta1 <- mcmc[,2]
  S <- nrow(mcmc)
  loocv <- numeric(n)
  for (i in 1:n){
    mu_i <- beta0 + beta1 * x[i]
    
    # p(y_i | theta_s)
    lik_i <- dnorm(y[i], mean = mu_i, sd = sqrt(sigma2))
    
    # Gelfand estimator
    loocv[i] <- 1 / mean(1 / lik_i)
  }
  
  elpd_rep[r] <- sum(log(loocv))
}

mean(elpd_rep)
var(elpd_rep)
