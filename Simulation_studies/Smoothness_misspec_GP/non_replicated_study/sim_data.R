# Data generation: y = X beta + w + eps, w ~ GP(0, sigma2 K_{phi,nu}), eps ~ N(0, tau2 I),
# r = sigma2 / (sigma2 + tau2); 20 x 20 regular grid on the unit square, X = (1, lat, lon).

# Matern correlation rho(h) = 2^{1-nu} / Gamma(nu) (h / phi)^nu K_nu(h / phi)
matern_cor <- function(D, phi, nu) {
  fields::Matern(D, range = phi, smoothness = nu)
}

sim_data <- function(n.grids, beta, sigma2, r, phi, nu, K, seed) {
  set.seed(seed)

  grid <- seq(0, 1, length.out = n.grids)
  locs <- as.matrix(expand.grid(lon = grid, lat = grid))
  N    <- nrow(locs)
  X    <- cbind(1, lat = locs[, "lat"], lon = locs[, "lon"])

  tau2  <- sigma2 * (1 - r) / r
  Sigma <- sigma2 * matern_cor(fields::rdist(locs), phi, nu)
  w     <- drop(crossprod(chol(Sigma), rnorm(N)))
  y     <- drop(X %*% beta) + w + rnorm(N, sd = sqrt(tau2))

  fold <- sample(rep(seq_len(K), length.out = N))

  list(y = y, X = X, locs = locs, fold = fold, K = K,
       truth = list(beta = beta, sigma2 = sigma2, tau2 = tau2, r = r, phi = phi, nu = nu))
}

if (sys.nframe() == 0) {
  if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

  sim_dat <- sim_data(n.grids = 10, 
                      beta = rep(0, 3), 
                      sigma2 = 1, 
                      r = 0.8, 
                      phi = 0.25,
                      nu = 0.5, 
                      K = 10, 
                      seed = 35)

  save(sim_dat, file = "sim_data.RData")
  cat(sprintf("Saved sim_data.RData: N = %d, K = %d, true nu = %.1f, r = %.2f\n",
              length(sim_dat$y), sim_dat$K, sim_dat$truth$nu, sim_dat$truth$r))
}
