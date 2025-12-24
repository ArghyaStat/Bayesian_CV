rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

# source("loocv_comparison.R")

load("elpd_outputs.RData")


###### Outputs ######

#### SE from our method (IS / MCSE-based) ####

cat("\nGelfand ELPD:", elpd_outputs$elpd_gelfand_is, 
    " SE:", elpd_outputs$se_elpd_gelfand_is, "\n")

cat("Zanella ELPD:", elpd_outputs$elpd_zanella_is, 
    " SE:", elpd_outputs$se_elpd_zanella_is, "\n")

cat("PSIS ELPD:", elpd_outputs$elpd_psis_is, 
    " SE:", elpd_outputs$se_elpd_psis_is, "\n")


#### SE from Monte Carlo estimator (replications) ####

cat("\nReplication results:\n")

cat("\nGelfand:\n")
cat("Mean ELPD:", mean(elpd_outputs$elpd_gelfand_rep), "\n")
cat("Var  ELPD:", var(elpd_outputs$elpd_gelfand_rep), "\n")

cat("\nZanella:\n")
cat("Mean ELPD:", mean(elpd_outputs$elpd_zanella_rep), "\n")
cat("Var  ELPD:", var(elpd_outputs$elpd_zanella_rep), "\n")

cat("\nPSIS:\n")
cat("Mean ELPD:", mean(elpd_outputs$elpd_psis_rep), "\n")
cat("Var  ELPD:", var(elpd_outputs$elpd_psis_rep), "\n")