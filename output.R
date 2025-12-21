rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

# source("loocv_comparison.R")

load("elpd_outputs.RData")


###### Outputs ######

#### SE from our method ####
cat("\nGelfand ELPD:", elpd_outputs$elpd_gelfand_is, 
    " SE:", elpd_outputs$se_elpd_gelfand_is, "\n")

cat("Zanella ELPD:", elpd_outputs$elpd_zanella_is, 
    " SE:", elpd_outputs$se_elpd_zanella_is, "\n")


#### SE from Monte Carlo estimator ####
cat("\nReplication results:\n")

cat("Gelfand:\n")
cat("Mean ELPD:", mean(elpd_outputs$elpd_gelfand_rep), "\n")
cat("Var  ELPD:", var(elpd_outputs$elpd_gelfand_rep), "\n")

cat("Zanella:\n")
cat("Mean ELPD:", mean(elpd_outputs$elpd_zanella_rep), "\n")
cat("Var  ELPD:", var(elpd_outputs$elpd_zanella_rep), "\n")
