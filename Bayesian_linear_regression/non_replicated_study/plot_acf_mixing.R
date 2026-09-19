## ==============================================================================
## plot_acf_mixing.R
## Non-replicated Bayesian Linear Regression Study: n = 500, d = 250
## Autocorrelation Function (ACF) Comparison:

rm(list = ls())
gc()

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

# 1. Load Outputs
if (!file.exists("M1_outputs.RData")) stop("M1_outputs.RData not found.")
load("M1_outputs.RData")
m1_acf_beta <- acf_worst_beta
m1_acf_elpd <- acf_elppd
m1_beta_idx <- worst_beta_idx
m1_S        <- max(history$S)
m1_ess_beta <- history$minESS_xi[length(history$minESS_xi)]
m1_ess_elpd <- history$ess_elppd[length(history$ess_elppd)]
rm(history, acc_mix, elpd_M1, elpd_M2, elpd_diff_val, acf_worst_beta, acf_elppd)

has_m2 <- file.exists("M2_outputs.RData")
if (has_m2) {
  load("M2_outputs.RData")
  m2_acf_beta <- acf_worst_beta
  m2_acf_elpd <- acf_elppd
  m2_beta_idx <- worst_beta_idx
  m2_S        <- max(history$S)
  m2_ess_beta <- history$minESS_xi[length(history$minESS_xi)]
  m2_ess_elpd <- history$ess_elppd[length(history$ess_elppd)]
  rm(history, acc_mix, elpd_M1, elpd_M2, elpd_diff_val, acf_worst_beta, acf_elppd)
}

# 2. Styling Palette
COL_BETA  <- "#2166AC"   # deep royal blue for worst beta component
COL_ELPD  <- "#D73027"   # deep crimson red for ELPD_hat
COL_BAND  <- "grey70"   # Bartlett CI bounds

pdf_file <- "acf_comparison.pdf"
pdf(pdf_file, width = if (has_m2) 18 else 10, height = 13, useDingbats = FALSE)

if (has_m2) {
  par(
    mfrow    = c(2, 2),
    mar      = c(5.5, 6.8, 3.2, 1.8),
    oma      = c(0.4, 0.4, 0.4, 0.4),
    mgp      = c(4.0, 1.3, 0),
    cex.lab  = 2.0,
    cex.axis = 1.6,
    las      = 1
  )
} else {
  par(
    mfrow    = c(2, 1),
    mar      = c(5.5, 6.8, 3.2, 1.8),
    oma      = c(0.4, 0.4, 0.4, 0.4),
    mgp      = c(4.0, 1.3, 0),
    cex.lab  = 2.0,
    cex.axis = 1.6,
    las      = 1
  )
}

max_l <- 50
disp_lags <- 0:max_l
idx <- disp_lags + 1

# Helper function to plot Row 1: ACF
plot_single_acf <- function(acf_beta, acf_elpd, beta_idx, ess_beta, ess_elpd, S_draws, panel_title) {
  y_beta   <- acf_beta[idx]
  y_elpd   <- acf_elpd[idx]
  ci_bound <- 1.96 / sqrt(S_draws)
  
  plot(NA, xlim = c(0, max_l), ylim = c(-0.15, 1.05),
       xlab = expression("MCMC Lag (" * italic(k) * ")"),
       ylab = expression("Autocorrelation (" * italic(rho) * ")"),
       main = panel_title, cex.main = 1.8, xaxt = "n", yaxt = "n")
  
  axis(1, at = seq(0, max_l, by = 10))
  axis(2, at = seq(-0.2, 1.0, by = 0.2))
  
  abline(h = seq(-0.2, 1.0, by = 0.2), col = "grey92", lty = 1, lwd = 0.8)
  abline(v = seq(0, max_l, by = 10), col = "grey92", lty = 1, lwd = 0.8)
  
  # Zero reference and Bartlett 95% CI bands
  abline(h = 0, col = "grey40", lty = 1, lwd = 1.2)
  abline(h = c(-ci_bound, ci_bound), col = COL_BAND, lty = 2, lwd = 1.8)
  
  # Worst Beta component
  lines(disp_lags, y_beta, col = COL_BETA, lwd = 3.2)
  points(disp_lags, y_beta, col = COL_BETA, pch = 19, cex = 1.1)
  
  # ELPD_hat process
  lines(disp_lags, y_elpd, col = COL_ELPD, lwd = 3.2)
  points(disp_lags, y_elpd, col = COL_ELPD, pch = 15, cex = 1.1)
  
  legend("topright", inset = c(0.02, 0.02),
         legend = c(
           bquote(widehat(ELPPD) ~ "( " * xi * ", " * widehat(ESS) == .(round(ess_elpd)) * ")"),
           bquote(beta[.(beta_idx)] ~ "(worst, " * widehat(ESS) == .(round(ess_beta)) * ")")
         ),
         col = c(COL_ELPD, COL_BETA, COL_BAND),
         pch = c(15, 19, NA),
         lty = c(1, 1, 2),
         lwd = c(3.2, 3.2, 1.8),
         pt.cex = c(1.1, 1.1, NA),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.2, y.intersp = 1.2)
  
  box(lwd = 1.8)
}

# Helper function to plot Row 2: Cumulative tau(K) = 1 + 2 * sum(rho)
plot_single_tau <- function(acf_beta, acf_elpd, beta_idx, ess_beta, ess_elpd, S_draws, panel_title) {
  tau_beta <- 1 + 2 * cumsum(acf_beta[2:(max_l + 1)])
  tau_elpd <- 1 + 2 * cumsum(acf_elpd[2:(max_l + 1)])
  y_max    <- max(c(tau_beta, tau_elpd), na.rm = TRUE) * 1.18
  
  plot(NA, xlim = c(1, max_l), ylim = c(1, y_max),
       xlab = expression("Truncation lag (" * italic(K) * ")"),
       ylab = expression("Integrated autocorr time " * widehat(tau)(italic(K))),
       main = panel_title, cex.main = 1.7, xaxt = "n", yaxt = "n")
  
  axis(1, at = seq(0, max_l, by = 10))
  axis(2, at = seq(1, ceiling(y_max), by = 1))
  
  abline(h = seq(1, ceiling(y_max), by = 1), col = "grey92", lty = 1, lwd = 0.8)
  abline(v = seq(0, max_l, by = 10), col = "grey92", lty = 1, lwd = 0.8)
  
  # Lines for cumulative tau
  lines(1:max_l, tau_beta, col = COL_BETA, lwd = 3.2)
  points(1:max_l, tau_beta, col = COL_BETA, pch = 19, cex = 1.0)
  
  lines(1:max_l, tau_elpd, col = COL_ELPD, lwd = 3.2)
  points(1:max_l, tau_elpd, col = COL_ELPD, pch = 15, cex = 1.0)
  
  legend("topleft", inset = c(0.02, 0.02),
         legend = c(
           bquote(widehat(tau)[ELPPD](italic(K))),
           bquote(widehat(tau)[beta](italic(K)))
         ),
         col = c(COL_ELPD, COL_BETA),
         pch = c(15, 19),
         lty = c(1, 1),
         lwd = c(3.2, 3.2),
         pt.cex = c(1.0, 1.0),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.2, y.intersp = 1.2)
  
  box(lwd = 1.8)
}

# Row 1, Column 1: M1 Autocorrelation
plot_single_acf(m1_acf_beta, m1_acf_elpd, m1_beta_idx, m1_ess_beta, m1_ess_elpd, m1_S,
                expression(bold("(a) Model ") * bolditalic(M)[1] ~ bold("(True): Autocorrelation ") * italic(rho)(italic(k))))

# Row 1, Column 2: M2 Autocorrelation (if available)
if (has_m2) {
  plot_single_acf(m2_acf_beta, m2_acf_elpd, m2_beta_idx, m2_ess_beta, m2_ess_elpd, m2_S,
                  expression(bold("(b) Model ") * bolditalic(M)[2] ~ bold("(Misspecified): Autocorrelation ") * italic(rho)(italic(k))))
}

# Row 2, Column 1: M1 Cumulative tau(K)
plot_single_tau(m1_acf_beta, m1_acf_elpd, m1_beta_idx, m1_ess_beta, m1_ess_elpd, m1_S,
                expression(bold("(c) Model ") * bolditalic(M)[1] ~ bold(": Cumulative ") * widehat(tau)(italic(K)) == 1 + 2 * sum(italic(rho)(italic(k)), italic(k)==1, italic(K))))

# Row 2, Column 2: M2 Cumulative tau(K) (if available)
if (has_m2) {
  plot_single_tau(m2_acf_beta, m2_acf_elpd, m2_beta_idx, m2_ess_beta, m2_ess_elpd, m2_S,
                  expression(bold("(d) Model ") * bolditalic(M)[2] ~ bold(": Cumulative ") * widehat(tau)(italic(K)) == 1 + 2 * sum(italic(rho)(italic(k)), italic(k)==1, italic(K))))
}

dev.off()
cat(sprintf("Saved enhanced 2-row ACF & tau comparison plot to %s\n", pdf_file))
