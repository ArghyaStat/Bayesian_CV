# Exact 10-fold ELJPD reference for M1 (nu = 0.5) and M2 (nu = 1) under the Vecchia model.
# With (beta, s2) integrated out, p(y_I | phi, r) is closed form for any index set I, so
#   p(y_{D_k} | y_{-D_k}) = p(y) / p(y_{-D_k}),   p(y_I) = int p(y_I | phi, r) p(phi) p(r) dphi dr,
# computed by quadrature on (log phi, logit r). Two grid sizes check convergence.

rm(list = ls())
invisible(gc())

libraries <- c("fields", "FNN", "parallel")
invisible(lapply(libraries, library, character.only = TRUE))

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

source("vecchia.R")
source("model_utils.R")
source("sim_data.R")

models  <- c(M1 = 0.5, M2 = 1)
m       <- 20
grids   <- list(coarse = c(180, 140), fine = c(360, 280))   # (n_phi, n_r)
n_cores <- max(1L, min(20L, detectCores() - 1L))

inp <- load_inputs("sim_data.RData", m = m)
N   <- length(inp$y)
K   <- inp$K
D   <- rdist(inp$locs)
nbs <- lapply(seq_len(N), function(i) setdiff(na.omit(inp$NNarray[i, ]), i))
sets <- c(list(seq_len(N)), lapply(seq_len(K), function(k) which(inp$fold != k)))

log_sum_exp <- function(x) max(x) + log(sum(exp(x - max(x))))

# Vecchia correlation (U'U)^{-1} for C = r R_phi + (1 - r) I
vecchia_cor <- function(phi, r, nu) {
  C <- r * matern_cor(D, phi, nu)
  diag(C) <- 1
  U <- matrix(0, N, N)
  for (i in seq_len(N)) {
    nb <- nbs[[i]]
    d  <- 1
    if (length(nb)) {
      b <- solve(C[nb, nb], C[nb, i])
      d <- 1 - sum(C[i, nb] * b)
      U[i, nb] <- -b / sqrt(d)
    }
    U[i, i] <- 1 / sqrt(d)
  }
  tcrossprod(forwardsolve(U, diag(N)))
}

# log p(y_I | phi, r): y_I ~ N(0, s2 (C_I + X_I V X_I')), s2 ~ IG(a, b) integrated out
log_marg <- function(C, idx) {
  y  <- inp$y[idx]
  X  <- inp$X[idx, , drop = FALSE]
  n  <- length(idx)
  a  <- inp$a_sigma
  b  <- inp$b_sigma
  Rc <- chol(C[idx, idx] + X %*% inp$V_prior %*% t(X))
  yt <- backsolve(Rc, y, transpose = TRUE)
  lgamma(a + n / 2) - lgamma(a) + a * log(b) - n / 2 * log(2 * pi) -
    sum(log(diag(Rc))) - (a + n / 2) * log(b + sum(yt^2) / 2)
}

exact_eljpd <- function(nu, n_phi, n_r) {
  u <- seq(log(0.002), log(0.9999 * inp$b_phi), length.out = n_phi)
  v <- seq(qlogis(0.002), qlogis(0.9995), length.out = n_r)
  g <- expand.grid(u = u, v = v)
  phi <- exp(g$u)
  r   <- plogis(g$v)
  # uniform priors on (0, b_phi) and (0, 1), times the Jacobian of (u, v) -> (phi, r)
  log_w <- log(diff(u[1:2]) * diff(v[1:2])) + g$u + log(r * (1 - r)) - log(inp$b_phi)

  L <- simplify2array(mclapply(seq_len(nrow(g)), function(j) {
    C <- vecchia_cor(phi[j], r[j], nu)
    vapply(sets, log_marg, numeric(1), C = C)
  }, mc.cores = n_cores))

  log_ml <- apply(L + matrix(log_w, nrow(L), ncol(L), byrow = TRUE), 1, log_sum_exp)
  fold   <- log_ml[1] - log_ml[-1]
  list(eljpd = sum(fold), fold = fold)
}

res <- lapply(models, function(nu) lapply(grids, function(gr) exact_eljpd(nu, gr[1], gr[2])))

for (mod in names(models))
  cat(sprintf("%s (nu = %.1f): ELJPD = %.4f (coarse), %.4f (fine)\n", mod, models[[mod]],
              res[[mod]]$coarse$eljpd, res[[mod]]$fine$eljpd))

exact <- lapply(res, function(x) list(eljpd = x$fine$eljpd, fold = x$fine$fold,
                                      grid_error = abs(x$fine$eljpd - x$coarse$eljpd)))
cat(sprintf("Exact gap M1 - M2: %.4f\n", exact$M1$eljpd - exact$M2$eljpd))

save(exact, models, m, grids, file = "exact_reference.RData")
