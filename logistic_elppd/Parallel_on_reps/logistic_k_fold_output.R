
rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

# -----------------------------
# LOAD DATA (KEEP SEPARATE)
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
# COLORS (by fold)
# -----------------------------
cols <- c("steelblue", "darkorange", "forestgreen")

# -----------------------------
# PREPARE DATA LISTS
# -----------------------------
# Order: 5-fold, 10-fold, n-fold

G_list <- list(res5[,"var_g"], res10[,"var_g"], resn[,"var_g"])
B_list <- list(res5[,"var_b"], res10[,"var_b"], resn[,"var_b"])
Z_list <- list(res5[,"var_z"], res10[,"var_z"], resn[,"var_z"])

# replicate each 3 times (for spacing groups)
expand_list <- function(L){
  list(L[[1]], L[[2]], L[[3]])
}

# -----------------------------
# PDF OUTPUT
# -----------------------------
pdf("kfold_9boxplots.pdf", width = 14, height = 4)

par(mfrow = c(1,3),
    mar  = c(3.2, 4, 2, 1),
    mgp  = c(2,0.7,0),
    tcl  = -0.3,
    cex.lab  = 1.5,
    cex.axis = 1.3,
    cex.main = 1.4)

# -----------------------------
# FUNCTION TO DRAW PANEL
# -----------------------------
plot_panel <- function(data_list, title){
  
  # positions for grouped spacing
  pos <- c(1,2,3, 5,6,7, 9,10,11)
  
  boxplot(
    c(data_list, data_list, data_list), # placeholder (fixed below)
    plot = FALSE
  )
  
  # build actual list (repeat each fold 3 times)
  dat <- list(
    data_list[[1]], data_list[[1]], data_list[[1]],
    data_list[[2]], data_list[[2]], data_list[[2]],
    data_list[[3]], data_list[[3]], data_list[[3]]
  )
  
  boxplot(dat,
          at = pos,
          col = rep(cols, each = 3),
          border = "black",
          xaxt = "n",
          main = title,
          ylab = "Variance")
  
  axis(1,
       at = c(2,6,10),
       labels = c("K=5","K=10","K=n"))
  
  abline(v = c(4,8), lty=2, col="gray70")
}

# -----------------------------
# PLOTS
# -----------------------------
plot_panel(G_list, "G (Gelfand)")
plot_panel(B_list, "B (Bronze)")
plot_panel(Z_list, "Z (Zanella)")

# legend (once)
legend("topright",
       legend = c("5-fold","10-fold","n-fold"),
       fill = cols,
       bty = "n",
       cex = 1.1)

dev.off()