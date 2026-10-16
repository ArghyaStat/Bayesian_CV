# Model M1 (Joint Covariance): K-fold ELJPD via Gibbs sampling from q_mix

rm(list = ls())
invisible(gc())
set.seed(35)

libraries <- c("fields", "FNN", "coda", "mcmcse", "sp", "gstat", "this.path",
               "Rcpp", "RcppArmadillo")
invisible(lapply(libraries, library, character.only = TRUE))

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}


source("vecchia.R")
source("aux_functions_joint.R")
source("eljpd.utils.R")
Rcpp::sourceCpp("gibbs_kfold_joint.cpp")

# ---- Parameters --------------------------------------------------------------

k.folds    <- 10
m          <- 20
t.index    <- 140
S_star     <- 5e3
S_max      <- 1e5
S_step     <- 5e3
S_vals     <- seq(S_star, S_max, by = S_step)
tuning.phi <- 1e-2

# ---- Data & Model Specs ------------------------------------------------------

data.candidates <- c("data_full.RData", "../data_full.RData", "../../data_full.RData")
data.file <- data.candidates[file.exists(data.candidates)][1]
load(data.file, eva.data <- new.env())
full.data.DF <- eva.data$data_DF

N.s <- 3503
spatial.subset <- seq(N.s * (t.index - 1) + 1, N.s * t.index)
data <- full.data.DF[spatial.subset, ]

locations <- cbind(data$lon, data$lat)
X <- cbind(1, locations)
Y <- cbind(log(1 + data$BA), data$CNT)
N <- nrow(Y); p <- ncol(X); q <- ncol(Y)
family <- c("Gaussian", "Poisson")

fold_id <- sample(rep(seq_len(k.folds), length.out = N))

# Vecchia ordering across all locations
ord      <- max_min(locations)
locs.ord <- locations[ord, , drop = FALSE]
X.ord    <- X[ord, , drop = FALSE]
Y.ord    <- Y[ord, , drop = FALSE]
fold.ord <- fold_id[ord]

NNarray <- neighbor_matrix(locs.ord, m)
storage.mode(NNarray) <- "integer"

b_phi <- max(rdist(locs.ord)) / 3

M.prior  <- matrix(0, p, q)
V.prior  <- 1e2 * diag(p)
S.prior  <- diag(q)
df.prior <- q + 1

W.init <- cbind(Y.ord[, 1], log(1 + Y.ord[, 2]))
Sigma  <- cov(W.init)
beta   <- solve(crossprod(X.ord)) %*% crossprod(X.ord, W.init)

# Variogram-based initialization for phi and nu
R.chol   <- chol(Sigma)
W.white  <- W.init %*% solve(R.chol)
df.joint <- data.frame(lon = rep(locs.ord[, 1], q), lat = rep(locs.ord[, 2], q),
                       values = as.vector(W.white))
coordinates(df.joint) <- ~ lon + lat
emp.vario.joint <- variogram(values ~ 1, df.joint)
init.vgm <- vgm(model = "Mat", psill = var(df.joint$values),
                range = median(emp.vario.joint$dist), nugget = 0)
fit.vgm <- fit.variogram(emp.vario.joint, model = init.vgm, fit.sills = TRUE)
phi <- fit.vgm$range[2]
nu  <- fit.vgm$kappa[2]

cat(sprintf("Model M1 setup: N = %d, q = %d, K = %d, phi0 = %.4f, nu = %.3f\n",
            N, q, k.folds, phi, nu))

# ---- Gibbs Sampling on q_mix -------------------------------------------------

gibbs_args <- list(
  Y = Y.ord, X = X.ord, locs = locs.ord, NNarray = NNarray,
  fold = as.integer(fold.ord), K = k.folds,
  family = family, family_par = rep(NA_real_, q), nu = nu,
  M_prior = M.prior, V_prior = V.prior, S_prior = S.prior, df_prior = df.prior,
  b_phi = b_phi, tuning_phi = tuning.phi
)

cat(sprintf("Running Gibbs sampler on q_mix: S = %d...\n", S_max))

t_start <- Sys.time()
fit_main <- do.call(gibbs_joint, c(gibbs_args, list(
  phi_init = phi, beta_init = beta, Sigma_init = Sigma,
  W_init = W.init, niters = S_max
)))
t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

acc.phi   <- fit_main$acc
theta     <- fit_main$samples
fold_ll   <- fit_main$fold_loglik
cat(sprintf("Gibbs sampling complete in %.2f mins. Acceptance rate: %.4f\n",
            t_elapsed, acc.phi))
cat(sprintf("Distinct folds visited: %d of %d\n", length(unique(fit_main$folds)), k.folds))
rm(fit_main)
invisible(gc())

# ---- Sequential ELJPD History ------------------------------------------------

cat("Evaluating sequential ELJPD and MCSE checkpoints...\n")
history <- list(
  S         = numeric(),
  eljpd_hat = numeric(),
  asymp_se  = numeric(),
  ess_elpd  = numeric(),
  minESS_xi = numeric()
)

for (S_curr in S_vals) {
  tup <- tuple_function(fold_ll, S_curr = S_curr)
  out <- eljpd.details(
    tup$numerator,
    tup$denominator,
    log_c = tup$log_c,
    samples = theta[1:S_curr, apply(theta[1:S_curr, , drop = FALSE], 2, sd) > 0, drop = FALSE]
  )
  history$S         <- c(history$S,         S_curr)
  history$eljpd_hat <- c(history$eljpd_hat, out$eljpd_hat)
  history$asymp_se  <- c(history$asymp_se,  out$asymp_se)
  history$ess_elpd  <- c(history$ess_elpd,  out$ess_elpd)
  history$minESS_xi <- c(history$minESS_xi, out$minESS_xi)
}

print(data.frame(
  S         = history$S,
  eljpd_hat = round(history$eljpd_hat, 3),
  asymp_se  = round(history$asymp_se, 4),
  ess_elpd  = round(history$ess_elpd, 1),
  minESS_xi = round(history$minESS_xi, 1)
))

# ---- Save --------------------------------------------------------------------

save(
  history,
  acc.phi,
  k.folds,
  m,
  S_star,
  S_max,
  file = "M1_outputs.RData"
)
cat("Results saved to M1_outputs.RData\n")
