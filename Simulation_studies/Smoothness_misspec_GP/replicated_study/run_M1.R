# Replicated study, model M1 (nu = 0.5): 500 independent q_mix chains on fixed data

rm(list = ls())
invisible(gc())

libraries <- c("fields", "FNN", "mcmcse", "foreach", "doParallel", "Rcpp", "RcppArmadillo")
invisible(lapply(libraries, library, character.only = TRUE))

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

source("vecchia.R")
source("eljpd.utils.R")
source("model_utils.R")
Rcpp::sourceCpp("gibbs_kfold.cpp")

model      <- "M1"
nu         <- 0.5
m          <- 20
tuning.phi <- 20
tuning.r   <- 20
S_star <- 5e3
S_max  <- 1.5e5
S_step <- 5e3
S_vals     <- seq(S_star, S_max, by = S_step)
reps       <- 500
seed0      <- 1904
n_cores    <- 50

inp <- load_inputs("sim_data.RData", m = m)

# FORK workers inherit the compiled sampler; it is not exported (its native
# pointer cannot be serialised)
cl <- tryCatch({
  c <- makeCluster(n_cores, type = "FORK")
  registerDoParallel(c)
  c
}, error = function(e) {
  registerDoSEQ()
  NULL
})

t_start <- Sys.time()
results_reps <- foreach(r = 1:reps, .packages = c("mcmcse"),
                        .noexport = c("gibbs_kfold", "matern")) %dopar% {

  set.seed(seed0 + r)
  fit <- run_gibbs(inp, nu = nu, niters = S_max, tuning_phi = tuning.phi,
                   tuning_r = tuning.r)

  list(history    = eljpd_history(fit, S_vals),
       acc_phi    = fit$acc_phi,
       acc_r      = fit$acc_r,
       fold_visit = tabulate(fit$folds, nbins = inp$K) / S_max)
}
t_mins <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

if (!is.null(cl)) stopCluster(cl)
cat(sprintf("%s: %d replicates in %.1f mins\n", model, reps, t_mins))

save(results_reps, model, nu, m, tuning.phi, tuning.r, S_vals, reps, seed0,
     t_mins, file = sprintf("%s_outputs_reps.RData", model))
