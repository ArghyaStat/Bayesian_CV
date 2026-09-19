## ==============================================================================
## output_nonreps.R
## Bayesian Linear Regression Case Study: n = 500, d = 250 (Mixture Gibbs Sampler)
##
## Non-replicated Study Analysis:
##   * Dynamically learns S_star, S_max, S_step from M1_outputs.RData & M2_outputs.RData
##   * Data-adaptive x-axis and y-axis limits/ticks for ELPPD and ESS
##   * Figure 1 (Side-by-Side):
##       - Panel (a): Dynamic ELPPD tracking with pointwise 95% asymptotic CIs,
##                    exact analytical reference lines, and empirical separation point S_sep
##       - Panel (b): Single-chain Effective Sample Size (ESS) trajectories for M1 & M2
##                    relative to sequential stopping boundaries (eps = 0.05 and 0.025)
##   * Figure 1 Standalone: High-resolution single-panel ELPPD tracking & separation
##   * Figure 2 (Efficiency Tracking):
##       - Panel (a): Raw effective sample sizes ESS_elppd and minESS_beta
##       - Panel (b): Relative efficiency ratios widehat(ESS) / S
##   * High-resolution PNG previews generated alongside publication PDFs
## ==============================================================================

rm(list = ls())
gc()

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

# ==============================================================================
# 1. Load Non-Replicated Datasets and Learn Grid Parameters
# ==============================================================================

if (!file.exists("M1_outputs.RData")) stop("M1_outputs.RData not found in current directory.")
if (!file.exists("M2_outputs.RData")) stop("M2_outputs.RData not found in current directory.")

load("M1_outputs.RData")
nr_M1_history <- history
nr_elpd_M1    <- elpd_M1
nr_elpd_M2    <- elpd_M2
acc_M1        <- acc_mix
d_val         <- d
n_val         <- n
eta_val       <- eta
rm(history, acc_mix, elpd_diff_val)

load("M2_outputs.RData")
nr_M2_history <- history
acc_M2        <- acc_mix
rm(history, acc_mix, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# Dynamically learn S_star, S_max, S_step from the data:
S_vals <- nr_M1_history$S
S_star <- min(S_vals)
S_max  <- max(S_vals)
S_step <- if (length(S_vals) > 1) diff(S_vals)[1] else NA_real_
n_checkpoints <- length(S_vals)

cat("==================================================================\n")
cat("Non-Replicated Output Analysis (Mixture Gibbs Sampler)\n")
cat("==================================================================\n")
cat(sprintf("Dimension d = %d, Sample size n = %d, Signal eta = %.2f\n", d_val, n_val, eta_val))
cat(sprintf("Iterations: S_star = %s, S_max = %s, S_step = %s (%d checkpoints)\n",
            format(S_star, big.mark = ","),
            format(S_max,  big.mark = ","),
            format(S_step, big.mark = ","),
            n_checkpoints))
cat(sprintf("Exact LOO-ELPD : M1 = %.4f, M2 = %.4f (Delta = %.4f)\n",
            nr_elpd_M1, nr_elpd_M2, nr_elpd_M1 - nr_elpd_M2))
cat(sprintf("HMC acceptance : M1 = %.4f, M2 = %.4f\n\n", acc_M1, acc_M2))

# ==============================================================================
# 2. Asymptotic Confidence Intervals & Data-Adaptive Separation Detection
# ==============================================================================

z_crit        <- qnorm(0.975)             # 1.959964
ess_target05  <- 4 * z_crit^2 / (0.05^2)   # ~6146.33
ess_target025 <- 4 * z_crit^2 / (0.025^2)  # ~24585.34

h1 <- nr_M1_history
h2 <- nr_M2_history

hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))

m1_lo <- h1$elppd_hat - hw1
m1_hi <- h1$elppd_hat + hw1
m2_lo <- h2$elppd_hat - hw2
m2_hi <- h2$elppd_hat + hw2

## Data-adaptive separation point S_sep:
## Earliest checkpoint S after which lower bound of M1 > upper bound of M2 through S_max
compute_S_sep <- function(h1, h2, z_crit = qnorm(0.975)) {
  hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
  hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))
  m1_lo <- h1$elppd_hat - hw1
  m2_hi <- h2$elppd_hat + hw2
  is_sep <- (m1_lo > m2_hi)
  if (length(is_sep) > 0 && is_sep[length(is_sep)]) {
    fi <- which(!is_sep)
    k_sep <- if (length(fi) == 0) 1L else max(fi) + 1L
    return(h1$S[k_sep])
  }
  return(NA_real_)
}

sep_S <- compute_S_sep(h1, h2, z_crit = z_crit)
if (!is.na(sep_S)) {
  cat(sprintf("Data-adaptive separation attained at S_sep = %s iterations!\n",
              format(sep_S, big.mark = ",")))
  k_sep <- which(h1$S == sep_S)
  cat(sprintf("  At S = %s: M1 lower bound = %.2f > M2 upper bound = %.2f (gap = +%.2f)\n\n",
              format(sep_S, big.mark = ","), m1_lo[k_sep], m2_hi[k_sep],
              m1_lo[k_sep] - m2_hi[k_sep]))
} else {
  cat("No complete separation attained through S_max.\n\n")
}

# ==============================================================================
# 3. Data-Adaptive Axis Limits, Ticks, and Colour Palettes
# ==============================================================================

COL_M1_A    <- "#08519C"   # deep royal blue (M1 true)
COL_M2_A    <- "#D73027"   # deep crimson red (M2 misspecified)
COL_SEP     <- "#006837"   # dark forest green (separation line)

COL_ESS_ELPD <- "#D73027"  # deep crimson red (ESS_elppd)
COL_MIN_ESS  <- "#2166AC"  # deep royal blue (minESS_beta)

# Data-adaptive x-axis limits and tick marks:
xlim_common      <- range(S_vals)
s_ticks_uniform  <- pretty(S_vals, n = 8)
s_ticks_uniform  <- s_ticks_uniform[s_ticks_uniform >= S_star & s_ticks_uniform <= S_max]
s_labels_uniform <- paste0(s_ticks_uniform / 1000, "K")

# Data-adaptive y-axis limits and ticks for ELPPD:
# Learned from range of ELPPD_hat and exact reference values across both models
elppd_min_raw <- min(c(h1$elppd_hat, h2$elppd_hat, nr_elpd_M1, nr_elpd_M2, m1_lo, m2_lo))
elppd_max_raw <- max(c(h1$elppd_hat, h2$elppd_hat, nr_elpd_M1, nr_elpd_M2))
y_range       <- elppd_max_raw - elppd_min_raw

# pretty ticks encompassing the data with top clearance for the legend:
y_ticks1   <- pretty(c(elppd_min_raw - 0.05 * y_range, elppd_max_raw + 0.35 * y_range), n = 7)
ylim_elppd <- range(y_ticks1)

# Data-adaptive y-axis limits and ticks for ESS:
ess_all_vals <- c(h1$ess_elppd, h1$minESS_xi, h2$ess_elppd, h2$minESS_xi, ess_target05, ess_target025)
max_ess_val  <- max(ess_all_vals, na.rm = TRUE)
y_ticks_ess  <- pretty(c(0, max_ess_val * 1.08), n = 5)
ylim_ess     <- c(0, max(y_ticks_ess))
y_labels_ess <- ifelse(y_ticks_ess == 0, "0",
                       ifelse(y_ticks_ess >= 1e3, paste0(y_ticks_ess / 1000, "K"), as.character(y_ticks_ess)))

# Data-adaptive y-axis limits and ticks for Relative Efficiency (ESS / S):
ratio_all_vals <- c(h1$minESS_xi / h1$S, h1$ess_elppd / h1$S, h2$minESS_xi / h2$S, h2$ess_elppd / h2$S)
max_ratio_val  <- max(ratio_all_vals, na.rm = TRUE)
y_ticks_ratio  <- pretty(c(0, max_ratio_val * 1.12), n = 6)
ylim_ratio     <- c(0, max(y_ticks_ratio))
y_labels_ratio <- sprintf("%.2f", y_ticks_ratio)

# ==============================================================================
# 4. FIGURE 1: 2-Panel Side-by-Side (a) ELPPD Tracking and (b) ESS Progression
# ==============================================================================

plot_side_by_side <- function() {
  par(
    mfrow    = c(1, 2),
    mar      = c(6.5, 8.5, 2.5, 2.0),
    oma      = c(0.5, 0.5, 0.5, 0.5),
    mgp      = c(5.2, 1.5, 0),
    cex.lab  = 2.5,
    cex.axis = 1.8,
    las      = 1
  )

  # ----------------------------------------------------------------------------
  # PANEL (a): Dynamic ELPPD Tracking & Point of Separation
  # ----------------------------------------------------------------------------
  plot(NA, xlim = xlim_common, ylim = ylim_elppd,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ELPPD)(italic(S))),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks1, labels = y_ticks1)

  # Subtle background grid
  abline(h = y_ticks1,        col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  # Exact LOO reference lines
  abline(h = nr_elpd_M1, col = COL_M1_A, lty = 2, lwd = 3.2)
  abline(h = nr_elpd_M2, col = COL_M2_A, lty = 2, lwd = 3.2)

  # Separation point vertical line
  if (!is.na(sep_S)) {
    abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.5)
  }

  # 95% Confidence interval arrows (top of M2 arrow clamped to ylim[2] to stay inside box)
  arrows(h1$S, m1_lo, h1$S, m1_hi, code = 3, angle = 90, length = 0.07,
         col = adjustcolor(COL_M1_A, 0.50), lwd = 2.6)
  arrows(h2$S, m2_lo, h2$S, pmin(m2_hi, ylim_elppd[2]), code = 3, angle = 90, length = 0.07,
         col = adjustcolor(COL_M2_A, 0.55), lwd = 2.6)

  # Trajectories: M1 solid filled circle (pch=19), M2 open circle (pch=1)
  lines(h1$S,  h1$elppd_hat, col = COL_M1_A, lwd = 3.5)
  points(h1$S, h1$elppd_hat, col = COL_M1_A, pch = 19, cex = 1.8)

  lines(h2$S,  h2$elppd_hat, col = COL_M2_A, lwd = 3.5)
  points(h2$S, h2$elppd_hat, col = COL_M2_A, pch = 1,  cex = 1.8, lwd = 2.2)

  # Separation label
  sep_legend_txt <- if (!is.na(sep_S)) {
    bquote("Separation (" * italic(S)[sep] == .(format(sep_S, big.mark = ",")) * ")")
  } else {
    expression("No separation")
  }

  legend("topright", inset = c(0.025, 0.025),
         legend = c(
           expression(italic(M)[1] ~ "(True)"),
           expression(italic(M)[2] ~ "(Misspecified)"),
           expression(ELPPD(italic(M)[1]) ~ "(Exact)"),
           expression(ELPPD(italic(M)[2]) ~ "(Exact)"),
           sep_legend_txt
         ),
         col = c(COL_M1_A, COL_M2_A, COL_M1_A, COL_M2_A, COL_SEP),
         pch = c(19, 1, NA, NA, NA),
         lty = c(1, 1, 2, 2, 4),
         lwd = c(3.2, 3.2, 3.2, 3.2, 3.5),
         pt.cex = c(1.8, 1.8, NA, NA, NA),
         pt.lwd = c(NA, 2.2, NA, NA, NA),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.50, y.intersp = 1.15)

  box(lwd = 1.8)

  # ----------------------------------------------------------------------------
  # PANEL (b): Single-Chain Effective Sample Size (ESS) Progression
  # ----------------------------------------------------------------------------
  plot(NA, xlim = xlim_common, ylim = ylim_ess,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ESS)(italic(S))),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks_ess, labels = y_labels_ess)

  abline(h = y_ticks_ess,     col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  # Threshold lines for sequential stopping rules
  abline(h = ess_target05,  col = "darkorange2", lty = 3, lwd = 2.8)
  abline(h = ess_target025, col = "purple3",     lty = 3, lwd = 2.8)

  # Separation line vertical
  if (!is.na(sep_S)) {
    abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.2)
  }

  # M1 ESS curves
  lines(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, lwd = 3.2, lty = 1)
  points(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, pch = 19, cex = 1.6)

  lines(h1$S, h1$minESS_xi, col = COL_MIN_ESS, lwd = 3.2, lty = 1)
  points(h1$S, h1$minESS_xi, col = COL_MIN_ESS, pch = 19, cex = 1.6)

  # M2 ESS curves
  lines(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, lwd = 3.0, lty = 2)
  points(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, pch = 1, cex = 1.6, lwd = 2.0)

  lines(h2$S, h2$minESS_xi, col = COL_MIN_ESS, lwd = 3.0, lty = 2)
  points(h2$S, h2$minESS_xi, col = COL_MIN_ESS, pch = 1, cex = 1.6, lwd = 2.0)

  legend("topleft", inset = c(0.025, 0.025),
         legend = c(
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[1])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[1])),
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[2])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[2])),
           expression(ESS^"*" ~ (epsilon == 0.05 ~ ":" ~ 6 * "," * 146)),
           expression(ESS^"*" ~ (epsilon == 0.025 ~ ":" ~ 24 * "," * 585))
         ),
         col   = c(COL_ESS_ELPD, COL_MIN_ESS, COL_ESS_ELPD, COL_MIN_ESS, "darkorange2", "purple3"),
         pch   = c(19, 19, 1, 1, NA, NA),
         lty   = c(1, 1, 2, 2, 3, 3),
         lwd   = c(3.2, 3.2, 3.0, 3.0, 2.8, 2.8),
         pt.cex = c(1.6, 1.6, 1.6, 1.6, NA, NA),
         pt.lwd = c(NA, NA, 2.0, 2.0, NA, NA),
         bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex   = 1.35, y.intersp = 1.15)

  box(lwd = 1.8)
}

# Generate PDF & PNG for Side-by-Side Plot (elppd_side_by_side_d250_n500 & elppd_mcse)
for (base_name in c("elppd_side_by_side_d250_n500", "elppd_mcse")) {
  pdf(paste0(base_name, ".pdf"), width = 20, height = 9.5, useDingbats = FALSE)
  plot_side_by_side()
  dev.off()

  png(paste0(base_name, ".png"), width = 3000, height = 1425, res = 150)
  plot_side_by_side()
  dev.off()
  cat(sprintf("Saved %s.pdf and %s.png\n", base_name, base_name))
}

# ==============================================================================
# 5. FIGURE 1 STANDALONE: Single-Panel ELPPD Tracking & Separation Point
# ==============================================================================

plot_elppd_standalone <- function() {
  par(
    mar      = c(6.5, 8.5, 2.5, 2.0),
    oma      = c(0.5, 0.5, 0.5, 0.5),
    mgp      = c(5.2, 1.5, 0),
    cex.lab  = 2.4,
    cex.axis = 1.9,
    las      = 1
  )

  plot(NA, xlim = xlim_common, ylim = ylim_elppd,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ELPPD)(italic(S))),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks1, labels = y_ticks1)

  abline(h = y_ticks1,        col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  abline(h = nr_elpd_M1, col = COL_M1_A, lty = 2, lwd = 3.2)
  abline(h = nr_elpd_M2, col = COL_M2_A, lty = 2, lwd = 3.2)

  if (!is.na(sep_S)) {
    abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.5)
  }

  arrows(h1$S, m1_lo, h1$S, m1_hi, code = 3, angle = 90, length = 0.07,
         col = adjustcolor(COL_M1_A, 0.50), lwd = 2.6)
  arrows(h2$S, m2_lo, h2$S, pmin(m2_hi, ylim_elppd[2]), code = 3, angle = 90, length = 0.07,
         col = adjustcolor(COL_M2_A, 0.55), lwd = 2.6)

  lines(h1$S,  h1$elppd_hat, col = COL_M1_A, lwd = 3.5)
  points(h1$S, h1$elppd_hat, col = COL_M1_A, pch = 19, cex = 1.8)

  lines(h2$S,  h2$elppd_hat, col = COL_M2_A, lwd = 3.5)
  points(h2$S, h2$elppd_hat, col = COL_M2_A, pch = 1,  cex = 1.8, lwd = 2.2)

  sep_legend_txt <- if (!is.na(sep_S)) {
    bquote("Separation (" * italic(S)[sep] == .(format(sep_S, big.mark = ",")) * ")")
  } else {
    expression("No separation")
  }

  legend("topright", inset = c(0.025, 0.025),
         legend = c(
           expression(italic(M)[1] ~ "(True)"),
           expression(italic(M)[2] ~ "(Misspecified)"),
           expression(ELPPD(italic(M)[1]) ~ "(Exact)"),
           expression(ELPPD(italic(M)[2]) ~ "(Exact)"),
           sep_legend_txt
         ),
         col = c(COL_M1_A, COL_M2_A, COL_M1_A, COL_M2_A, COL_SEP),
         pch = c(19, 1, NA, NA, NA),
         lty = c(1, 1, 2, 2, 4),
         lwd = c(3.2, 3.2, 3.2, 3.2, 3.5),
         pt.cex = c(1.8, 1.8, NA, NA, NA),
         pt.lwd = c(NA, 2.2, NA, NA, NA),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.40, y.intersp = 1.15)

  box(lwd = 1.8)
}

pdf("elppd_tracking_d250_n500.pdf", width = 11, height = 9, useDingbats = FALSE)
plot_elppd_standalone()
dev.off()

png("elppd_tracking_d250_n500.png", width = 1650, height = 1350, res = 150)
plot_elppd_standalone()
dev.off()
cat("Saved elppd_tracking_d250_n500.pdf and elppd_tracking_d250_n500.png\n")

# ==============================================================================
# 6. FIGURE 2: MCMC Efficiency - (a) Raw ESS and (b) Relative Efficiency ESS / S
# ==============================================================================

plot_ess_2panel <- function() {
  par(
    mfrow    = c(1, 2),
    mar      = c(6.5, 8.5, 2.5, 2.0),
    oma      = c(0.5, 0.5, 0.5, 0.5),
    mgp      = c(5.2, 1.5, 0),
    cex.lab  = 2.5,
    cex.axis = 1.8,
    las      = 1
  )

  # Panel (a): Raw ESS
  plot(NA, xlim = xlim_common, ylim = ylim_ess,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ESS)),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks_ess, labels = y_labels_ess)

  abline(h = y_ticks_ess,     col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  # Stopping thresholds
  abline(h = ess_target05,  col = "darkorange2", lty = 3, lwd = 2.8)
  abline(h = ess_target025, col = "purple3",     lty = 3, lwd = 2.8)

  # Curves
  lines(h1$S, h1$minESS_xi, col = COL_MIN_ESS, lwd = 3.2, lty = 1)
  points(h1$S, h1$minESS_xi, col = COL_MIN_ESS, pch = 19, cex = 1.7)

  lines(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, lwd = 3.2, lty = 1)
  points(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, pch = 19, cex = 1.7)

  lines(h2$S, h2$minESS_xi, col = COL_MIN_ESS, lwd = 3.0, lty = 2)
  points(h2$S, h2$minESS_xi, col = COL_MIN_ESS, pch = 1, cex = 1.7, lwd = 2.0)

  lines(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, lwd = 3.0, lty = 2)
  points(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, pch = 1, cex = 1.7, lwd = 2.0)

  legend("topleft", inset = c(0.025, 0.025),
         legend = c(
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[1])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[1])),
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[2])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[2])),
           expression(ESS^"*" ~ (epsilon == 0.05 ~ ":" ~ 6 * "," * 146)),
           expression(ESS^"*" ~ (epsilon == 0.025 ~ ":" ~ 24 * "," * 585))
         ),
         col   = c(COL_ESS_ELPD, COL_MIN_ESS, COL_ESS_ELPD, COL_MIN_ESS, "darkorange2", "purple3"),
         pch   = c(19, 19, 1, 1, NA, NA),
         lty   = c(1, 1, 2, 2, 3, 3),
         lwd   = c(3.2, 3.2, 3.0, 3.0, 2.8, 2.8),
         pt.cex = c(1.7, 1.7, 1.7, 1.7, NA, NA),
         pt.lwd = c(NA, NA, 2.0, 2.0, NA, NA),
         bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex   = 1.35, y.intersp = 1.15)
  box(lwd = 1.8)

  # Panel (b): Relative Efficiency ESS / S
  plot(NA, xlim = xlim_common, ylim = ylim_ratio,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ESS) / italic(S)),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks_ratio, labels = y_labels_ratio)

  abline(h = y_ticks_ratio,   col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  ratio_m1_min <- h1$minESS_xi / h1$S
  ratio_m1_elp <- h1$ess_elppd / h1$S
  ratio_m2_min <- h2$minESS_xi / h2$S
  ratio_m2_elp <- h2$ess_elppd / h2$S

  lines(h1$S, ratio_m1_min, col = COL_MIN_ESS, lwd = 3.2, lty = 1)
  points(h1$S, ratio_m1_min, col = COL_MIN_ESS, pch = 19, cex = 1.7)

  lines(h1$S, ratio_m1_elp, col = COL_ESS_ELPD, lwd = 3.2, lty = 1)
  points(h1$S, ratio_m1_elp, col = COL_ESS_ELPD, pch = 19, cex = 1.7)

  lines(h2$S, ratio_m2_min, col = COL_MIN_ESS, lwd = 3.0, lty = 2)
  points(h2$S, ratio_m2_min, col = COL_MIN_ESS, pch = 1, cex = 1.7, lwd = 2.0)

  lines(h2$S, ratio_m2_elp, col = COL_ESS_ELPD, lwd = 3.0, lty = 2)
  points(h2$S, ratio_m2_elp, col = COL_ESS_ELPD, pch = 1, cex = 1.7, lwd = 2.0)

  legend("topright", inset = c(0.025, 0.025),
         legend = c(
           expression(widehat(ESS)[ELPPD] / italic(S) ~ (italic(M)[1])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) / italic(S) ~ (italic(M)[1])),
           expression(widehat(ESS)[ELPPD] / italic(S) ~ (italic(M)[2])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) / italic(S) ~ (italic(M)[2]))
         ),
         col   = c(COL_ESS_ELPD, COL_MIN_ESS, COL_ESS_ELPD, COL_MIN_ESS),
         pch   = c(19, 19, 1, 1),
         lty   = c(1, 1, 2, 2),
         lwd   = c(3.2, 3.2, 3.0, 3.0),
         pt.cex = c(1.7, 1.7, 1.7, 1.7),
         pt.lwd = c(NA, NA, 2.0, 2.0),
         bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex   = 1.35, y.intersp = 1.15)
  box(lwd = 1.8)
}

for (base_name in c("ess_tracking", "ess_boxplot")) {
  pdf(paste0(base_name, ".pdf"), width = 20, height = 9.5, useDingbats = FALSE)
  plot_ess_2panel()
  dev.off()

  png(paste0(base_name, ".png"), width = 3000, height = 1425, res = 150)
  plot_ess_2panel()
  dev.off()
  cat(sprintf("Saved %s.pdf and %s.png\n", base_name, base_name))
}

# Standalone Raw ESS
plot_ess_raw_standalone <- function() {
  par(
    mar      = c(6.2, 7.8, 1.5, 1.5),
    oma      = c(0.2, 0.2, 0.2, 0.2),
    mgp      = c(4.8, 1.4, 0),
    cex.lab  = 2.2,
    cex.axis = 1.7,
    las      = 1
  )

  plot(NA, xlim = xlim_common, ylim = ylim_ess,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ESS)),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks_ess, labels = y_labels_ess)

  abline(h = y_ticks_ess,     col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  abline(h = ess_target05,  col = "darkorange2", lty = 3, lwd = 2.5)
  abline(h = ess_target025, col = "purple3",     lty = 3, lwd = 2.5)

  lines(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, lwd = 3.0, lty = 1)
  points(h1$S, h1$ess_elppd, col = COL_ESS_ELPD, pch = 19, cex = 1.5)

  lines(h1$S, h1$minESS_xi, col = COL_MIN_ESS, lwd = 3.0, lty = 1)
  points(h1$S, h1$minESS_xi, col = COL_MIN_ESS, pch = 19, cex = 1.5)

  lines(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, lwd = 2.8, lty = 2)
  points(h2$S, h2$ess_elppd, col = COL_ESS_ELPD, pch = 1, cex = 1.5, lwd = 1.8)

  lines(h2$S, h2$minESS_xi, col = COL_MIN_ESS, lwd = 2.8, lty = 2)
  points(h2$S, h2$minESS_xi, col = COL_MIN_ESS, pch = 1, cex = 1.5, lwd = 1.8)

  legend("topleft", inset = c(0.04, 0.04),
         legend = c(
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[1])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[1])),
           expression(widehat(ESS)[ELPPD] ~ (italic(M)[2])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) ~ (italic(M)[2])),
           expression(ESS^"*" ~ (epsilon == 0.05)),
           expression(ESS^"*" ~ (epsilon == 0.025))
         ),
         col   = c(COL_ESS_ELPD, COL_MIN_ESS, COL_ESS_ELPD, COL_MIN_ESS, "darkorange2", "purple3"),
         pch   = c(19, 19, 1, 1, NA, NA),
         lty   = c(1, 1, 2, 2, 3, 3),
         lwd   = c(3.0, 3.0, 2.8, 2.8, 2.5, 2.5),
         pt.cex = c(1.5, 1.5, 1.5, 1.5, NA, NA),
         bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex   = 1.25, y.intersp = 1.15)
  box(lwd = 1.8)
}

for (base_name in c("ess_raw", "ess_raw_boxplot")) {
  pdf(paste0(base_name, ".pdf"), width = 10, height = 8, useDingbats = FALSE)
  plot_ess_raw_standalone()
  dev.off()

  png(paste0(base_name, ".png"), width = 1500, height = 1200, res = 150)
  plot_ess_raw_standalone()
  dev.off()
}
cat("Saved ess_raw.pdf / ess_raw_boxplot.pdf\n")

# Standalone Relative Efficiency
plot_ess_ratio_standalone <- function() {
  par(
    mar      = c(6.2, 7.8, 1.5, 1.5),
    oma      = c(0.2, 0.2, 0.2, 0.2),
    mgp      = c(4.8, 1.4, 0),
    cex.lab  = 2.2,
    cex.axis = 1.7,
    las      = 1
  )

  plot(NA, xlim = xlim_common, ylim = ylim_ratio,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = expression(widehat(ESS) / italic(S)),
       main = "", xaxt = "n", yaxt = "n")

  axis(1, at = s_ticks_uniform, labels = s_labels_uniform, gap.axis = 0)
  axis(2, at = y_ticks_ratio, labels = y_labels_ratio)

  abline(h = y_ticks_ratio,   col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks_uniform, col = "grey92", lty = 1, lwd = 0.8)

  ratio_m1_min <- h1$minESS_xi / h1$S
  ratio_m1_elp <- h1$ess_elppd / h1$S
  ratio_m2_min <- h2$minESS_xi / h2$S
  ratio_m2_elp <- h2$ess_elppd / h2$S

  lines(h1$S, ratio_m1_elp, col = COL_ESS_ELPD, lwd = 3.0, lty = 1)
  points(h1$S, ratio_m1_elp, col = COL_ESS_ELPD, pch = 19, cex = 1.5)

  lines(h1$S, ratio_m1_min, col = COL_MIN_ESS, lwd = 3.0, lty = 1)
  points(h1$S, ratio_m1_min, col = COL_MIN_ESS, pch = 19, cex = 1.5)

  lines(h2$S, ratio_m2_elp, col = COL_ESS_ELPD, lwd = 2.8, lty = 2)
  points(h2$S, ratio_m2_elp, col = COL_ESS_ELPD, pch = 1, cex = 1.5, lwd = 1.8)

  lines(h2$S, ratio_m2_min, col = COL_MIN_ESS, lwd = 2.8, lty = 2)
  points(h2$S, ratio_m2_min, col = COL_MIN_ESS, pch = 1, cex = 1.5, lwd = 1.8)

  legend("topright", inset = c(0.04, 0.04),
         legend = c(
           expression(widehat(ESS)[ELPPD] / italic(S) ~ (italic(M)[1])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) / italic(S) ~ (italic(M)[1])),
           expression(widehat(ESS)[ELPPD] / italic(S) ~ (italic(M)[2])),
           expression(widehat(ESS) ~ "of" ~ min(beta[italic(j)]) / italic(S) ~ (italic(M)[2]))
         ),
         col   = c(COL_ESS_ELPD, COL_MIN_ESS, COL_ESS_ELPD, COL_MIN_ESS),
         pch   = c(19, 19, 1, 1),
         lty   = c(1, 1, 2, 2),
         lwd   = c(3.0, 3.0, 2.8, 2.8),
         pt.cex = c(1.5, 1.5, 1.5, 1.5),
         bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex   = 1.25, y.intersp = 1.15)
  box(lwd = 1.8)
}

for (base_name in c("ess_ratio", "ess_ratio_boxplot")) {
  pdf(paste0(base_name, ".pdf"), width = 10, height = 8, useDingbats = FALSE)
  plot_ess_ratio_standalone()
  dev.off()

  png(paste0(base_name, ".png"), width = 1500, height = 1200, res = 150)
  plot_ess_ratio_standalone()
  dev.off()
}
cat("Saved ess_ratio.pdf / ess_ratio_boxplot.pdf\n")

# ==============================================================================
# 7. Print Checkpoint Summary Table
# ==============================================================================

cat("\n==================================================================\n")
cat("Selected Checkpoints Summary Table\n")
cat("==================================================================\n")
sel_idx <- sort(unique(c(1L,
                         which(h1$S %in% c(20e3, 50e3, 85e3, 100e3, 150e3, 200e3)),
                         n_checkpoints)))

summary_df <- data.frame(
  S         = format(h1$S[sel_idx], big.mark = ","),
  M1_ELPPD  = sprintf("%.2f [%.2f, %.2f]", h1$elppd_hat[sel_idx], m1_lo[sel_idx], m1_hi[sel_idx]),
  M2_ELPPD  = sprintf("%.2f [%.2f, %.2f]", h2$elppd_hat[sel_idx], m2_lo[sel_idx], m2_hi[sel_idx]),
  Separated = ifelse(m1_lo[sel_idx] > m2_hi[sel_idx], "YES (M1>M2)", "No"),
  M1_ESS_EL = format(round(h1$ess_elppd[sel_idx]), big.mark = ","),
  M1_minESS = format(round(h1$minESS_xi[sel_idx]), big.mark = ","),
  M2_ESS_EL = format(round(h2$ess_elppd[sel_idx]), big.mark = ","),
  M2_minESS = format(round(h2$minESS_xi[sel_idx]), big.mark = ",")
)

print(summary_df, row.names = FALSE)
cat("==================================================================\n")
cat("Analysis complete.\n")
