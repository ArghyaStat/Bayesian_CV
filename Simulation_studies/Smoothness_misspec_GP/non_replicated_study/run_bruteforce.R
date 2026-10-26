# Brute-force 10-fold ELJPD reference for M1 and M2: one chain on p(theta | y_{-D_k}) per fold k

rm(list = ls())
invisible(gc())

libraries <- c("fields", "FNN", "mcmcse", "parallel", "Rcpp", "RcppArmadillo")
invisible(lapply(libraries, library, character.only = TRUE))

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

source("vecchia.R")
source("eljpd.utils.R")
source("model_utils.R")
Rcpp::sourceCpp("gibbs_kfold.cpp")

models  <- c(M1 = 0.5, M2 = 1)
tuning  <- list(M1 = c(1.3, 5.5), M2 = c(0.5, 2.0))
m       <- 20
S_bf    <- 5e4
n_cores <- max(1L, min(20L, detectCores() - 1L))

inp  <- load_inputs("sim_data.RData", m = m)
jobs <- expand.grid(k = seq_len(inp$K), model = names(models), stringsAsFactors = FALSE)

fold_fits <- mclapply(seq_len(nrow(jobs)), function(j) {
  set.seed(5000 + j)
  k   <- jobs$k[j]
  fit <- run_gibbs(inp, nu = models[[jobs$model[j]]], niters = S_bf,
                   tuning_phi = tuning[[jobs$model[j]]][1],
                   tuning_r = tuning[[jobs$model[j]]][2], fixed_fold = k)
  ll  <- fit$fold_loglik[, k]
  w   <- exp(ll - max(ll))
  c(elpd_k = max(ll) + log(mean(w)),
    se_k   = mcse(w, method = "bm")$se / mean(w))
}, mc.cores = n_cores)
fold_fits <- do.call(rbind, fold_fits)

bf <- lapply(names(models), function(mod) {
  idx <- jobs$model == mod
  list(fold_elpd = fold_fits[idx, "elpd_k"],
       eljpd     = sum(fold_fits[idx, "elpd_k"]),
       se        = sqrt(sum(fold_fits[idx, "se_k"]^2)))
})
names(bf) <- names(models)

cat(sprintf("Brute-force ELJPD  M1: %.3f (%.3f)  M2: %.3f (%.3f)  diff: %.3f\n",
            bf$M1$eljpd, bf$M1$se, bf$M2$eljpd, bf$M2$se, bf$M1$eljpd - bf$M2$eljpd))

save(bf, S_bf, tuning, m, file = "bruteforce_reference.RData")
