## ==============================================================================
## output_analysis.R
## Bayesian Linear Regression Case Study: n = 500, d = 250 (Mixture IS)
##
## Panel (a): Single-chain ELPPD tracking with 95% CI and separation point S_sep
##            using updated non-replicated outputs up to S_max = 250,000.
## Panel (b): Replicated termination distributions (S_term, ELPPD_hat) for
##            eps = 0.05 & 0.025 across 500 independent MCMC runs.
##            - solid round points with distinct harmonious shades for M1 and M2
##            - eps = 0.05: deep crimson red (M1) vs soft coral-rose red (M2)
##            - eps = 0.025: deep royal blue (M1) vs soft sky blue (M2)
## Figure 2 : Superimposed boxplots across S iterations:
##            - Panel (a): Raw effective sample sizes ESS_elpd and minESS_xi
##            - Panel (b): Relative efficiency ratios widehat(ESS)/S
## ==============================================================================

rm(list = ls())
gc()

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

S_star <- 2e4
S_max  <- 2.5e5
S_step <- 5e3
S_vals <- seq(S_star, S_max, by = S_step)  # 47 checkpoints
reps   <- 500

z_crit        <- qnorm(0.975)             # 1.959964
ess_target05  <- 4 * z_crit^2 / (0.05^2)   # ~6146.33
ess_target025 <- 4 * z_crit^2 / (0.025^2)  # ~24585.34


load("M1_outputs.RData")
nr_M1_history <- history
nr_elpd_M1    <- elpd_M1
nr_elpd_M2    <- elpd_M2
rm(history, acc_mix, elpd_diff_val)

load("M2_outputs.RData")
nr_M2_history <- history
rm(history, acc_mix, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# Replicated runs (500 reps, S_star = 20K to S_max = 250K)
load("M1_outputs_reps.RData")
rep_M1      <- results_reps
rep_elpd_M1 <- elpd_M1
rep_elpd_M2 <- elpd_M2
rm(results_reps, elpd_diff_val)

load("M2_outputs_reps.RData")
rep_M2 <- results_reps
rm(results_reps, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# 2. Helper Functions: Termination and Data-Adaptive Separation Detection

get_term <- function(h, ess_target) {
  idx_ok <- which(h$ess_elppd >= ess_target)
  k      <- if (length(idx_ok) > 0) idx_ok[1] else length(h$S)
  list(S_term    = h$S[k],
       ELPPD_hat = h$elppd_hat[k],
       asymp_se  = h$asymp_se[k],
       k         = k)
}

## Data-adaptive separation point S_sep for single chain tracking:
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

## Data-adaptive separation point across replicated chains:
## Supports median of per-replicate S_sep or checkpoint-wise empirical quantile envelopes
compute_replicated_S_sep <- function(rep_M1, rep_M2, z_crit = qnorm(0.975),
                                     method = c("median_per_rep", "quantile_envelope"),
                                     alpha = 0.05) {
  method <- match.arg(method)
  if (method == "median_per_rep") {
    sep_per_rep <- sapply(seq_along(rep_M1), function(i) {
      compute_S_sep(rep_M1[[i]]$history, rep_M2[[i]]$history, z_crit = z_crit)
    })
    return(list(
      S_sep   = median(sep_per_rep, na.rm = TRUE),
      mean    = mean(sep_per_rep, na.rm = TRUE),
      summary = summary(sep_per_rep),
      all     = sep_per_rep
    ))
  } else if (method == "quantile_envelope") {
    S_vals <- rep_M1[[1]]$history$S
    n_S <- length(S_vals)
    q_lo <- rep(NA_real_, n_S)
    q_hi <- rep(NA_real_, n_S)
    for (i in seq_along(S_vals)) {
      s <- S_vals[i]
      v1 <- sapply(rep_M1, function(r) {
        idx <- which(r$history$S == s)
        if (length(idx)) r$history$elppd_hat[idx] else NA_real_
      })
      v2 <- sapply(rep_M2, function(r) {
        idx <- which(r$history$S == s)
        if (length(idx)) r$history$elppd_hat[idx] else NA_real_
      })
      q_lo[i] <- quantile(v1, alpha / 2, na.rm = TRUE)
      q_hi[i] <- quantile(v2, 1 - alpha / 2, na.rm = TRUE)
    }
    is_sep <- (q_lo > q_hi)
    sep_S <- NA_real_
    if (length(is_sep) > 0 && is_sep[length(is_sep)]) {
      fi <- which(!is_sep)
      k_sep <- if (length(fi) == 0) 1L else max(fi) + 1L
      sep_S <- S_vals[k_sep]
    }
    return(list(S_sep = sep_S, q_lo = q_lo, q_hi = q_hi, S_vals = S_vals))
  }
}

t05_M1 <- lapply(rep_M1, function(x) get_term(x$history, ess_target05))
t025_M1 <- lapply(rep_M1, function(x) get_term(x$history, ess_target025))
t05_M2 <- lapply(rep_M2, function(x) get_term(x$history, ess_target05))
t025_M2 <- lapply(rep_M2, function(x) get_term(x$history, ess_target025))

t05_s_M1 <- sapply(t05_M1, `[[`, "S_term")
t05_y_M1 <- sapply(t05_M1, `[[`, "ELPPD_hat")
t025_s_M1 <- sapply(t025_M1, `[[`, "S_term")
t025_y_M1 <- sapply(t025_M1, `[[`, "ELPPD_hat")

t05_s_M2 <- sapply(t05_M2, `[[`, "S_term")
t05_y_M2 <- sapply(t05_M2, `[[`, "ELPPD_hat")
t025_s_M2 <- sapply(t025_M2, `[[`, "S_term")
t025_y_M2 <- sapply(t025_M2, `[[`, "ELPPD_hat")

# 3. Colour Palette Definitions
# ==============================================================================
# Panel (a): Deep royal blue (M1) and crimson red (M2)
COL_M1_A <- "#08519C"   # deep royal blue
COL_M2_A <- "#D73027"   # deep crimson red
COL_SEP  <- "#006837"   # forest green for separation

# Panel (b): Distinct solid shades (no open points, no wireframe clutter)
# eps = 0.05: deep crimson red for M1, rich golden-yellow for M2
COL_05_M1 <- "#D73027"   # deep crimson red
COL_05_M2 <- "#D4AC0D"   # rich golden-yellow

# eps = 0.025: deep royal blue for M1, soft sky blue for M2
COL_025_M1 <- "#08519C"   # deep royal blue
COL_025_M2 <- "#6BAED6"   # soft sky-blue tint

# ==============================================================================
# 4. FIGURE 1: Side-by-Side (a) Single Chain and (b) Replicated Scatter
# ==============================================================================
# Uniform x-axis ticks across all plots: 20K to 240K spaced by 20K
s_ticks_uniform  <- seq(20e3, 240e3, by = 20e3)
s_labels_uniform <- paste0(s_ticks_uniform / 1000, "K")

pdf_file1 <- "elppd_mcse.pdf"
pdf(pdf_file1, width = 20, height = 9.5, useDingbats = FALSE)

par(
  mfrow    = c(1, 2),
  mar      = c(6.5, 8.5, 2.5, 2.0),
  oma      = c(0.5, 0.5, 0.5, 0.5),
  mgp      = c(5.2, 1.5, 0),
  cex.lab  = 2.6,
  cex.axis = 2.1,
  las      = 1
)

# ------------------------------------------------------------------------------
# PANEL (a): Dynamic ELPPD Tracking & Point of Separation (S up to 250,000)
# ------------------------------------------------------------------------------
h1 <- nr_M1_history
h2 <- nr_M2_history

hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))

m1_lo <- h1$elppd_hat - hw1
m1_hi <- h1$elppd_hat + hw1
m2_lo <- h2$elppd_hat - hw2
m2_hi <- h2$elppd_hat + hw2

# Separation point: data-adaptive detection via compute_S_sep()
sep_S <- compute_S_sep(h1, h2, z_crit = z_crit)

ylim1 <- c(-855.0, -825.0)
xlim1 <- range(h1$S)

plot(NA, xlim = xlim1, ylim = ylim1,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ELPPD)(italic(S))),
     main = "", xaxt = "n", yaxt = "n")

s_ticks1  <- s_ticks_uniform
s_labels1 <- s_labels_uniform
axis(1, at = s_ticks1, labels = s_labels1)

y_ticks1 <- seq(-855, -825, by = 5)
axis(2, at = y_ticks1, labels = y_ticks1)

# Subtle background grid
abline(h = y_ticks1, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks1, col = "grey92", lty = 1, lwd = 0.8)

# Exact LOO reference lines
abline(h = nr_elpd_M1, col = COL_M1_A, lty = 2, lwd = 3.2)
abline(h = nr_elpd_M2, col = COL_M2_A, lty = 2, lwd = 3.2)

# Separation point vertical line
if (!is.na(sep_S)) {
  abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.5)
}

# 95% Confidence interval arrows
arrows(h1$S, m1_lo, h1$S, m1_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M1_A, 0.50), lwd = 2.8)
arrows(h2$S, m2_lo, h2$S, m2_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M2_A, 0.55), lwd = 2.8)

# Trajectories: M1 solid filled circle (pch=19), M2 open circle (pch=1)
lines(h1$S,  h1$elppd_hat, col = COL_M1_A, lwd = 3.5)
points(h1$S, h1$elppd_hat, col = COL_M1_A, pch = 19, cex = 1.8)

lines(h2$S,  h2$elppd_hat, col = COL_M2_A, lwd = 3.5)
points(h2$S, h2$elppd_hat, col = COL_M2_A, pch = 1,  cex = 1.8, lwd = 2.2)

# Legend placed in top-right
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
       cex = 1.55, y.intersp = 1.15)

box(lwd = 1.8)

# ------------------------------------------------------------------------------
# PANEL (b): Replicated Termination Distributions (S_term, ELPPD_hat)
# ------------------------------------------------------------------------------
xlim2 <- c(20e3, 250e3)
ylim2 <- c(-855.0, -825.0)

plot(NA, xlim = xlim2, ylim = ylim2,
     xlab = expression("Iterations for termination (" * italic(S)[term] * ")"),
     ylab = expression(widehat(ELPPD)(italic(S)[term])),
     main = "", xaxt = "n", yaxt = "n")

s_ticks2  <- s_ticks_uniform
s_labels2 <- s_labels_uniform
axis(1, at = s_ticks2, labels = s_labels2)

y_ticks2 <- seq(-855, -825, by = 5)
axis(2, at = y_ticks2, labels = y_ticks2)

# Subtle background grid
abline(h = y_ticks2, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks2, col = "grey92", lty = 1, lwd = 0.8)

# Exact LOO reference lines
abline(h = rep_elpd_M1, col = "grey35", lty = 2, lwd = 3.0)
abline(h = rep_elpd_M2, col = "grey35", lty = 2, lwd = 3.0)

# eps = 0.05: Crimson (M1) and Golden-yellow (M2)
points(t05_s_M2, t05_y_M2,
       col = adjustcolor(COL_05_M2, 0.75), pch = 19, cex = 1.05)
points(t05_s_M1, t05_y_M1,
       col = adjustcolor(COL_05_M1, 0.70), pch = 19, cex = 1.05)

# eps = 0.025: Blue family
points(t025_s_M2, t025_y_M2,
       col = adjustcolor(COL_025_M2, 0.65), pch = 19, cex = 1.05)
points(t025_s_M1, t025_y_M1,
       col = adjustcolor(COL_025_M1, 0.70), pch = 19, cex = 1.05)

# Legend placed cleanly at bottom-right (termination distributions)
legend("topright", inset = c(0.025, 0.025),
       legend = c(
         expression(epsilon == 0.05 * ":" ~ italic(M)[1] ~ "(True)"),
         expression(epsilon == 0.05 * ":" ~ italic(M)[2] ~ "(Misspecified)"),
         expression(epsilon == 0.025 * ":" ~ italic(M)[1] ~ "(True)"),
         expression(epsilon == 0.025 * ":" ~ italic(M)[2] ~ "(Misspecified)"),
         expression("Exact" ~ ELPPD ~ "(Reference)")
       ),
       col   = c(COL_05_M1, COL_05_M2, COL_025_M1, COL_025_M2, "grey35"),
       pch   = c(19, 19, 19, 19, NA),
       lty   = c(NA, NA, NA, NA, 2),
       lwd   = c(NA, NA, NA, NA, 3.0),
       pt.cex = c(1.5, 1.5, 1.5, 1.5, NA),
       bty   = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex   = 1.35, y.intersp = 1.15)

## Note on data-adaptive replicated S_sep:
## If one wishes to compute the separation point across the 500 replicated runs,
## call compute_replicated_S_sep(rep_M1, rep_M2, z_crit = z_crit):
##   - method = 'median_per_rep': median S_sep across runs = 60,000 (IQR: [35K, 88.8K])
##   - method = 'quantile_envelope': 95% envelope separation S_sep = 85,000
## Panel (b) displays the bivariate termination scatter (S_term, ELPPD_hat), which
## cleanly separates M1 and M2 in ELPPD without an arbitrary partition across S_term.

box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file1))

# 5. FIGURE 2: Superimposed Boxplots of ESS and Relative Efficiency (ESS / S)
# ==============================================================================

n_steps <- length(S_vals)
ess_elpd_mat   <- matrix(NA_real_, reps, n_steps)
min_ess_xi_mat <- matrix(NA_real_, reps, n_steps)

for (r in seq_len(reps)) {
  h <- rep_M1[[r]]$history
  m <- match(S_vals, h$S)
  valid <- !is.na(m)
  ess_elpd_mat[r, valid]   <- h$ess_elppd[m[valid]]
  min_ess_xi_mat[r, valid] <- h$minESS_xi[m[valid]]
}

# Compute relative efficiency ratios: ESS / S (current value of S)
ratio_elpd   <- sweep(ess_elpd_mat, 2, S_vals, "/")
ratio_min_xi <- sweep(min_ess_xi_mat, 2, S_vals, "/")

COL_ESS_ELPD <- "#D73027"   # deep crimson red
COL_MIN_ESS  <- "#2166AC"   # deep royal blue

box_w   <- S_step * 0.45
xlim_bp <- c(15e3, 255e3)

# Helper to draw boxplot
draw_bxp <- function(vals, at_pos, bw, fill_col, border_col) {
  vals <- vals[is.finite(vals)]
  if (length(vals) < 5L) return(invisible(NULL))
  bp <- boxplot(vals, plot = FALSE)
  bxp(bp, at = at_pos, add = TRUE, boxwex = bw,
      col = fill_col, border = border_col,
      outline = FALSE, axes = FALSE,
      whisklty = 1, staplelty = 1, medlwd = 2.4)
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 5A. 2-PANEL FIGURE: Side-by-Side (a) Raw ESS and (b) Relative Efficiency ESS/S
# ------------------------------------------------------------------------------
pdf_file2 <- "ess_boxplot.pdf"
pdf(pdf_file2, width = 20, height = 9.5, useDingbats = FALSE)

par(
  mfrow    = c(1, 2),
  mar      = c(6.5, 8.5, 2.5, 2.0),
  oma      = c(0.5, 0.5, 0.5, 0.5),
  mgp      = c(5.2, 1.5, 0),
  cex.lab  = 2.6,
  cex.axis = 2.1,
  las      = 1
)

# PANEL (a): Superimposed Raw ESS
ylim_bp <- c(0, 6e4)
plot(NA, xlim = xlim_bp, ylim = ylim_bp,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ESS)),
     main = "", xaxt = "n", yaxt = "n")

s_ticks3  <- s_ticks_uniform
s_labels3 <- s_labels_uniform
axis(1, at = s_ticks3, labels = s_labels3)
y_ticks3  <- seq(0, 50000, by = 10000)
y_labels3 <- ifelse(y_ticks3 == 0, "0", paste0(y_ticks3 / 1000, "K"))
axis(2, at = y_ticks3, labels = y_labels3)

abline(h = y_ticks3, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks3, col = "grey92", lty = 1, lwd = 0.8)

for (k in seq_len(n_steps)) {
  draw_bxp(min_ess_xi_mat[, k], S_vals[k], box_w,
           adjustcolor(COL_MIN_ESS,  0.40), COL_MIN_ESS)
  draw_bxp(ess_elpd_mat[, k], S_vals[k], box_w,
           adjustcolor(COL_ESS_ELPD, 0.40), COL_ESS_ELPD)
}

legend("topleft", inset = c(0.025, 0.025),
       legend = c(
         expression(min[italic(j)] ~ widehat(ESS)(theta[italic(j)])),
         expression(widehat(ESS)[ELPPD])
       ),
       fill   = c(adjustcolor(COL_MIN_ESS, 0.40), adjustcolor(COL_ESS_ELPD, 0.40)),
       border = c(COL_MIN_ESS, COL_ESS_ELPD),
       bty    = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex    = 1.45, y.intersp = 1.2)
box(lwd = 1.8)

# PANEL (b): Superimposed Relative Efficiency: ESS / S
ylim_ratio <- c(0, 0.38)
plot(NA, xlim = xlim_bp, ylim = ylim_ratio,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ESS) / italic(S)),
     main = "", xaxt = "n", yaxt = "n")

axis(1, at = s_ticks3, labels = s_labels3)
y_ticks_ratio  <- seq(0, 0.35, by = 0.05)
y_labels_ratio <- sprintf("%.2f", y_ticks_ratio)
axis(2, at = y_ticks_ratio, labels = y_labels_ratio)

abline(h = y_ticks_ratio, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks3,       col = "grey92", lty = 1, lwd = 0.8)

for (k in seq_len(n_steps)) {
  draw_bxp(ratio_min_xi[, k], S_vals[k], box_w,
           adjustcolor(COL_MIN_ESS,  0.40), COL_MIN_ESS)
  draw_bxp(ratio_elpd[, k], S_vals[k], box_w,
           adjustcolor(COL_ESS_ELPD, 0.40), COL_ESS_ELPD)
}

legend("topright", inset = c(0.025, 0.025),
       legend = c(
         expression(min[italic(j)] ~ widehat(ESS)(theta[italic(j)]) / italic(S)),
         expression(widehat(ESS)[ELPPD] / italic(S))
       ),
       fill   = c(adjustcolor(COL_MIN_ESS, 0.40), adjustcolor(COL_ESS_ELPD, 0.40)),
       border = c(COL_MIN_ESS, COL_ESS_ELPD),
       bty    = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex    = 1.45, y.intersp = 1.2)
box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file2))

# ------------------------------------------------------------------------------
# 5B. STANDALONE SINGLE-PANEL FIGURES (for individual inclusion in LaTeX)
# ------------------------------------------------------------------------------
# Standalone Raw ESS Plot
pdf_file2_raw <- "ess_raw_boxplot.pdf"
pdf(pdf_file2_raw, width = 10, height = 8, useDingbats = FALSE)
par(
  mar      = c(6.2, 7.8, 1.5, 1.5),
  oma      = c(0.2, 0.2, 0.2, 0.2),
  mgp      = c(4.8, 1.4, 0),
  cex.lab  = 2.2,
  cex.axis = 1.7,
  las      = 1
)
plot(NA, xlim = xlim_bp, ylim = ylim_bp,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ESS)),
     main = "", xaxt = "n", yaxt = "n")
axis(1, at = s_ticks3, labels = s_labels3)
axis(2, at = y_ticks3, labels = y_labels3)
abline(h = y_ticks3, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks3, col = "grey92", lty = 1, lwd = 0.8)
for (k in seq_len(n_steps)) {
  draw_bxp(min_ess_xi_mat[, k], S_vals[k], box_w,
           adjustcolor(COL_MIN_ESS,  0.40), COL_MIN_ESS)
  draw_bxp(ess_elpd_mat[, k], S_vals[k], box_w,
           adjustcolor(COL_ESS_ELPD, 0.40), COL_ESS_ELPD)
}
legend("topleft", inset = c(0.04, 0.04),
       legend = c(
         expression(min[italic(j)] ~ widehat(ESS)(theta[italic(j)])),
         expression(widehat(ESS)[ELPPD])
       ),
       fill   = c(adjustcolor(COL_MIN_ESS, 0.40), adjustcolor(COL_ESS_ELPD, 0.40)),
       border = c(COL_MIN_ESS, COL_ESS_ELPD),
       bty    = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex    = 1.35, y.intersp = 1.2)
box(lwd = 1.8)
dev.off()
cat(sprintf("Saved %s\n", pdf_file2_raw))

# Standalone Relative Efficiency Plot (ESS / S)
pdf_file2_ratio <- "ess_ratio_boxplot.pdf"
pdf(pdf_file2_ratio, width = 10, height = 8, useDingbats = FALSE)
par(
  mar      = c(6.2, 7.8, 1.5, 1.5),
  oma      = c(0.2, 0.2, 0.2, 0.2),
  mgp      = c(4.8, 1.4, 0),
  cex.lab  = 2.2,
  cex.axis = 1.7,
  las      = 1
)
plot(NA, xlim = xlim_bp, ylim = ylim_ratio,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ESS) / italic(S)),
     main = "", xaxt = "n", yaxt = "n")
axis(1, at = s_ticks3, labels = s_labels3)
axis(2, at = y_ticks_ratio, labels = y_labels_ratio)
abline(h = y_ticks_ratio, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks3,       col = "grey92", lty = 1, lwd = 0.8)
for (k in seq_len(n_steps)) {
  draw_bxp(ratio_min_xi[, k], S_vals[k], box_w,
           adjustcolor(COL_MIN_ESS,  0.40), COL_MIN_ESS)
  draw_bxp(ratio_elpd[, k], S_vals[k], box_w,
           adjustcolor(COL_ESS_ELPD, 0.40), COL_ESS_ELPD)
}
legend("topright", inset = c(0.04, 0.04),
       legend = c(
         expression(min[italic(j)] ~ widehat(ESS)(theta[italic(j)]) / italic(S)),
         expression(widehat(ESS)[ELPPD] / italic(S))
       ),
       fill   = c(adjustcolor(COL_MIN_ESS, 0.40), adjustcolor(COL_ESS_ELPD, 0.40)),
       border = c(COL_MIN_ESS, COL_ESS_ELPD),
       bty    = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex    = 1.35, y.intersp = 1.2)
box(lwd = 1.8)
dev.off()
cat(sprintf("Saved %s\n", pdf_file2_ratio))
