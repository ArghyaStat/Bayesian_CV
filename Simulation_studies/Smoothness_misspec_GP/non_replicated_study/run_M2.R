# Non-replicated study, model M2 (nu = 1): 10-fold ELJPD from a single q_mix chain

rm(list = ls())
invisible(gc())

libraries <- c("fields", "FNN", "mcmcse", "Rcpp", "RcppArmadillo")
invisible(lapply(libraries, library, character.only = TRUE))

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

source("vecchia.R")
source("eljpd.utils.R")
source("model_utils.R")
Rcpp::sourceCpp("gibbs_kfold.cpp")

model  <- "M2"

nu     <- 1
m      <- 20
tuning.phi <- 16
tuning.r   <- 12
S_star <- 5e3
S_max  <- 1.5e5
S_step <- 5e3

S_vals <- seq(S_star, S_max, by = S_step)

inp <- load_inputs("sim_data.RData", m = m)

set.seed(4343)
t_start <- Sys.time()
fit     <- run_gibbs(inp, nu = nu, niters = S_max, tuning_phi = tuning.phi,
                     tuning_r = tuning.r)
t_mins  <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("%s: %.2f mins, acceptance phi %.3f, r %.3f\n", model, t_mins, fit$acc_phi, fit$acc_r))

history    <- eljpd_history(fit, S_vals)
acc_phi    <- fit$acc_phi
acc_r      <- fit$acc_r
fold_visit <- tabulate(fit$folds, nbins = inp$K) / S_max

print(data.frame(S = history$S, 
                 eljpd_hat = round(history$eljpd_hat, 3),
                 asymp_se = round(history$asymp_se, 4), 
                 ess_elpd = round(history$ess_elpd),
                 minESS_xi = round(history$minESS_xi)))

save(history, acc_phi, acc_r, fold_visit, model, nu, m, tuning.phi, tuning.r, S_vals, t_mins,
     file = sprintf("%s_outputs.RData", model))
