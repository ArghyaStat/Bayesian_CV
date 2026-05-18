
rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

load("elppd_exact_logistic.RData")

#===========================================================
# Fold settings
#===========================================================

K.vec <- c(5, 10, 50)

#===========================================================
# Colors
#===========================================================

submethod.cols <- c(
  "steelblue",   # IS
  "darkorange",  # TIS
  "forestgreen"  # PSIS
)

box.cols <- rep(
  submethod.cols,
  3
)

#===========================================================
# Helper functions
#===========================================================

plot_elpd_panel <- function(
    results,
    true_elpd,
    title){
  
  plot.data <- list(
    
    results$elppd_g_is,
    results$elppd_g_tis,
    results$elppd_g_psis,
    
    results$elppd_b_is,
    results$elppd_b_tis,
    results$elppd_b_psis,
    
    results$elppd_z_is,
    results$elppd_z_tis,
    results$elppd_z_psis
  )
  
  boxplot(
    plot.data,
    
    names = rep(
      c("IS", "TIS", "PSIS"),
      3
    ),
    
    col = box.cols,
    border = "black",
    
    las = 2,
    outline = FALSE,
    xaxt = "n",
    
    ylab = "ELPPD",
    main = title
  )
  
  axis(
    1,
    at = 1:9,
    labels = rep(
      c("IS", "TIS", "PSIS"),
      3
    ),
    las = 2
  )
  
  axis(
    1,
    at = c(2, 5, 8),
    labels = c(
      "Gelfand",
      "Bronze",
      "Zanella"
    ),
    tick = FALSE,
    line = 3
  )
  
  abline(
    h = true_elpd,
    col = "red",
    lwd = 2,
    lty = 2
  )
  
  legend(
    "topright",
    legend = c(
      "IS",
      "TIS",
      "PSIS",
      "Exact ELPPD"
    ),
    
    fill = c(
      submethod.cols,
      NA
    ),
    
    border = c(
      "black",
      "black",
      "black",
      NA
    ),
    
    lty = c(
      NA,
      NA,
      NA,
      2
    ),
    
    lwd = c(
      NA,
      NA,
      NA,
      2
    ),
    
    col = c(
      NA,
      NA,
      NA,
      "red"
    ),
    
    bty = "n",
    cex = 0.9
  )
}


plot_var_panel <- function(
    results,
    title){
  
  #---------------------------------------------------------
  # Estimated asymptotic variances
  #---------------------------------------------------------
  
  plot.data <- list(
    
    results$var_g_is,
    results$var_g_tis,
    results$var_g_psis,
    
    results$var_b_is,
    results$var_b_tis,
    results$var_b_psis,
    
    results$var_z_is,
    results$var_z_tis,
    results$var_z_psis
  )
  
  #---------------------------------------------------------
  # Sample variances
  #---------------------------------------------------------
  
  samp.var <- c(
    
    var(results$elppd_g_is),
    var(results$elppd_g_tis),
    var(results$elppd_g_psis),
    
    var(results$elppd_b_is),
    var(results$elppd_b_tis),
    var(results$elppd_b_psis),
    
    var(results$elppd_z_is),
    var(results$elppd_z_tis),
    var(results$elppd_z_psis)
  )
  
  ymax <- max(
    unlist(plot.data),
    samp.var
  )
  
  pos <- 1:9
  
  boxplot(
    plot.data,
    
    at = pos,
    
    names = rep(
      c("IS", "TIS", "PSIS"),
      3
    ),
    
    col = box.cols,
    border = "black",
    
    ylim = c(0, ymax),
    
    las = 2,
    outline = FALSE,
    xaxt = "n",
    
    ylab = "Variance",
    main = title
  )
  
  axis(
    1,
    at = pos,
    labels = rep(
      c("IS", "TIS", "PSIS"),
      3
    ),
    las = 2
  )
  
  axis(
    1,
    at = c(2, 5, 8),
    labels = c(
      "Gelfand",
      "Bronze",
      "Zanella"
    ),
    tick = FALSE,
    line = 3
  )
  
  # sample variance benchmark
  
  segments(
    pos - 0.25,
    samp.var,
    pos + 0.25,
    samp.var,
    lwd = 3,
    col = "red"
  )
  
  legend(
    "topright",
    legend = c(
      "IS",
      "TIS",
      "PSIS",
      "Sample Variance"
    ),
    
    fill = c(
      submethod.cols,
      NA
    ),
    
    border = c(
      "black",
      "black",
      "black",
      NA
    ),
    
    lwd = c(
      NA,
      NA,
      NA,
      3
    ),
    
    col = c(
      NA,
      NA,
      NA,
      "red"
    ),
    
    bty = "n",
    cex = 0.9
  )
}


#===========================================================
# ELPD comparison plot
#===========================================================

pdf(
  "logistic_k_fold_elpd_comparison.pdf",
  width = 16,
  height = 5
)

par(
  mfrow = c(1,3),
  mar = c(6,4,3,1),
  mgp = c(2,0.7,0),
  tcl = -0.3
)

for(K in K.vec){
  
  load(
    paste0(
      K,
      "_fold_elppd_logistic.RData"
    )
  )
  
  plot_elpd_panel(
    results = results,
    true_elpd = elppd_exact,
    title = paste0("K = ", K)
  )
}

dev.off()

#===========================================================
# Variance comparison plot
#===========================================================

pdf(
  "logistic_k_fold_variance_comparison.pdf",
  width = 16,
  height = 5
)

par(
  mfrow = c(1,3),
  mar = c(6,4,3,1),
  mgp = c(2,0.7,0),
  tcl = -0.3
)

for(K in K.vec){
  
  load(
    paste0(
      K,
      "_fold_elppd_logistic.RData"
    )
  )
  
  plot_var_panel(
    results = results,
    title = paste0("K = ", K)
  )
}

dev.off()