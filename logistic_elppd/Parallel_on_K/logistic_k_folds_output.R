rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

# -----------------------------
# LOAD DATA
# -----------------------------
load("5_fold_elppd_logistic.RData")
res5 <- results

load("10_fold_elppd_logistic.RData")
res10 <- results

load("n_fold_elppd_logistic.RData")
resn <- results

# -----------------------------
# COLUMN NAMES
# -----------------------------
cn <- c(
  "elpd_g","var_g","varS_g",
  "elpd_b","var_b","varS_b",
  "elpd_z","var_z","varS_z",
  "acc_g","acc_b","acc_z"
)

colnames(res5)  <- cn
colnames(res10) <- cn
colnames(resn)  <- cn

# -----------------------------
# COLORS
# -----------------------------
cols <- c("steelblue", "darkorange", "forestgreen")

# -----------------------------
# BUILD DATA (CLT vs Sample)
# -----------------------------
build_method_list <- function(res5, res10, resn, var_name, varS_name){
  
  list(
    res5[,var_name],  res5[,varS_name],
    res10[,var_name], res10[,varS_name],
    resn[,var_name],  resn[,varS_name]
  )
}

G_list <- build_method_list(res5, res10, resn, "var_g", "varS_g")
B_list <- build_method_list(res5, res10, resn, "var_b", "varS_b")
Z_list <- build_method_list(res5, res10, resn, "var_z", "varS_z")

# -----------------------------
# PDF
# -----------------------------
pdf("kfold_variance_comparison.pdf", width = 14, height = 4)

par(mfrow = c(1,3),
    mar  = c(3.5, 4, 2, 1),
    mgp  = c(2,0.7,0),
    tcl  = -0.3,
    cex.lab  = 1.4,
    cex.axis = 1.2,
    cex.main = 1.4)

# -----------------------------
# PLOT FUNCTION
# -----------------------------
plot_panel <- function(data_list, title){
  
  # positions (pairwise spacing)
  pos <- c(1,2, 4,5, 7,8)
  
  boxplot(data_list,
          at = pos,
          col = rep(cols, each = 2),
          border = "black",
          xaxt = "n",
          main = title,
          ylab = "Variance")
  
  axis(1,
       at = c(1.5,4.5,7.5),
       labels = c("K=5","K=10","K=n"))
  
  # vertical separators
  abline(v = c(3,6), lty=2, col="gray70")
}

# -----------------------------
# PLOTS
# -----------------------------
plot_panel(G_list, "G (Gelfand)")
plot_panel(B_list, "B (Bronze)")
plot_panel(Z_list, "Z (Zanella)")

# -----------------------------
# LEGEND
# -----------------------------
legend("topright",
       legend = c("CLT variance","Sample variance"),
       fill = c("gray70","gray30"),
       bty = "n",
       cex = 1.1)

dev.off()