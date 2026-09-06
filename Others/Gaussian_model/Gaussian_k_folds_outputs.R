rm(list = ls())

mydir <- this.path::here()
setwd(mydir)


load("5_fold_gaussian_exact_elpd.RData")

res5 <- results
true5 <- true_elpd

load("10_fold_gaussian_exact_elpd.RData")

res10 <- results
true10 <- true_elpd

load("50_fold_gaussian_exact_elpd.RData")

resn <- results
truen <- true_elpd


cn <- c(
  "elpd_g","var_g","acc_g",
  "elpd_b","var_b","acc_b",
  "elpd_z","var_z","acc_z"
)

colnames(res5)  <- cn
colnames(res10) <- cn
colnames(resn)  <- cn



varS5 <- c(
  var(res5[,"elpd_g"]),
  var(res5[,"elpd_b"]),
  var(res5[,"elpd_z"])
)

varS10 <- c(
  var(res10[,"elpd_g"]),
  var(res10[,"elpd_b"]),
  var(res10[,"elpd_z"])
)

varSn <- c(
  var(resn[,"elpd_g"]),
  var(resn[,"elpd_b"]),
  var(resn[,"elpd_z"])
)

# Variance dataframe
data5 <- list(
  res5[,"var_g"],
  res5[,"var_b"],
  res5[,"var_z"]
)

data10 <- list(
  res10[,"var_g"],
  res10[,"var_b"],
  res10[,"var_z"]
)

datan <- list(
  resn[,"var_g"],
  resn[,"var_b"],
  resn[,"var_z"]
)

# ELPD dataframe

data5_elpd <- list(
  res5[,"elpd_g"],
  res5[,"elpd_b"],
  res5[,"elpd_z"]
)

data10_elpd <- list(
  res10[,"elpd_g"],
  res10[,"elpd_b"],
  res10[,"elpd_z"]
)

datan_elpd <- list(
  resn[,"elpd_g"],
  resn[,"elpd_b"],
  resn[,"elpd_z"]
)



cols <- c(
  "steelblue",
  "darkorange",
  "forestgreen"
)

# Variance plot

pdf("Gaussian_k_fold_variance_comparison.pdf",
    width = 12,
    height = 4)

par(mfrow = c(1,3),
    mar  = c(3.5,4,2,1),
    mgp  = c(2,0.7,0),
    tcl  = -0.3,
    cex.lab  = 1.3,
    cex.axis = 1.1,
    cex.main = 1.3)

plot_panel <- function(data_list,
                       sample_var,
                       title){
  
  pos <- c(1,2,3)
  
  boxplot(data_list,
          at = pos,
          col = cols,
          border = "black",
          xaxt = "n",
          ylab = "Variance",
          main = title)
  
  axis(1,
       at = pos,
       labels = c(
         "Gelfand",
         "Bronze",
         "Zanella"
       ))
  
  # sample variance benchmark
  segments(pos - 0.25,
           sample_var,
           pos + 0.25,
           sample_var,
           lwd = 3,
           col = "red")
  
  legend("topright",
         legend = "Sample variance",
         lwd = 3,
         col = "red",
         bty = "n",
         cex = 1)
}

plot_panel(
  data5,
  varS5,
  "K = 5"
)

plot_panel(
  data10,
  varS10,
  "K = 10"
)

plot_panel(
  datan,
  varSn,
  "K = n"
)

dev.off()

# ELPD plot

pdf("Gaussian_k_fold_elpd_comparison.pdf",
    width = 12,
    height = 4)

par(mfrow = c(1,3),
    mar  = c(3.5,4,2,1),
    mgp  = c(2,0.7,0),
    tcl  = -0.3,
    cex.lab  = 1.3,
    cex.axis = 1.1,
    cex.main = 1.3)

plot_elpd_panel <- function(data_list,
                            true_elpd,
                            title){
  
  pos <- c(1,2,3)
  
  boxplot(data_list,
          at = pos,
          col = cols,
          border = "black",
          xaxt = "n",
          ylab = "ELPD",
          main = title)
  
  axis(1,
       at = pos,
       labels = c(
         "Gelfand",
         "Bronze",
         "Zanella"
       ))
  
  abline(h = true_elpd,
         col = "red",
         lwd = 3)
  
  legend("topright",
         legend = "True ELPD",
         lwd = 3,
         col = "red",
         bty = "n",
         cex = 1)
}

plot_elpd_panel(
  data5_elpd,
  true5,
  "K = 5"
)

plot_elpd_panel(
  data10_elpd,
  true10,
  "K = 10"
)

plot_elpd_panel(
  datan_elpd,
  truen,
  "K = n"
)

dev.off()