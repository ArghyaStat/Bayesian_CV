## Subset selection in Bayesian linear regression (Section 6.1): M1 (true) vs M2 (one predictor dropped)
## elpd_mcse_linreg.pdf: (a) single-chain ELPPD tracking with 95% MC intervals and S_sep;
##                       (b) replicated termination points (S_term, ELPPD_hat) for eps = 0.05, 0.025

rm(list = ls())
gc()

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

z_crit        <- qnorm(0.975)             # 1.959964
ess_target05  <- 4 * z_crit^2 / (0.05^2)   # ~6146.33
ess_target025 <- 4 * z_crit^2 / (0.025^2)  # ~24585.34

nonrep_dir <- "../non_replicated_study"
rep_dir    <- "../replicated_study"

# Single-chain runs (S up to 250K)
load(file.path(nonrep_dir, "M1_outputs.RData"))
nr_M1_history <- history
nr_elpd_M1    <- elpd_M1
nr_elpd_M2    <- elpd_M2
rm(history, acc_mix, elpd_diff_val)

load(file.path(nonrep_dir, "M2_outputs.RData"))
nr_M2_history <- history
rm(history, acc_mix, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# Replicated runs (500 reps, S from 20K to 250K)
load(file.path(rep_dir, "M1_outputs_reps.RData"))
rep_M1      <- results_reps
rep_elpd_M1 <- elpd_M1
rep_elpd_M2 <- elpd_M2
rm(results_reps, elpd_diff_val)

load(file.path(rep_dir, "M2_outputs_reps.RData"))
rep_M2 <- results_reps
rm(results_reps, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# First checkpoint with ESS_ELPPD >= target (S_max if never reached)
get_term <- function(h, ess_target) {
  idx_ok <- which(h$ess_elppd >= ess_target)
  k      <- if (length(idx_ok) > 0) idx_ok[1] else length(h$S)
  list(S_term    = h$S[k],
       ELPPD_hat = h$elppd_hat[k],
       asymp_se  = h$asymp_se[k],
       k         = k)
}

# Earliest checkpoint after which the lower 95% bound of M1 stays above the upper bound of M2
# (asymp_se stores the asymptotic variance, hence the square root)
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

# Separation point across replicates: median per-replicate S_sep or quantile envelopes
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

# Colours: panel (a) by model, panel (b) by eps and model
COL_M1_A <- "#08519C"   # deep royal blue
COL_M2_A <- "#D73027"   # deep crimson red
COL_SEP  <- "#006837"   # forest green for separation
COL_05_M1 <- "#D73027"   # deep crimson red
COL_05_M2 <- "#D4AC0D"   # rich golden-yellow
COL_025_M1 <- "#08519C"   # deep royal blue
COL_025_M2 <- "#6BAED6"   # soft sky-blue tint

# x-axis ticks: 20K to 240K by 20K
s_ticks_uniform  <- seq(20e3, 240e3, by = 20e3)
s_labels_uniform <- paste0(s_ticks_uniform / 1000, "K")

pdf_file1 <- "elpd_mcse_linreg.pdf"
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

# Panel (a): single-chain ELPPD tracking and separation point
h1 <- nr_M1_history
h2 <- nr_M2_history

hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))

m1_lo <- h1$elppd_hat - hw1
m1_hi <- h1$elppd_hat + hw1
m2_lo <- h2$elppd_hat - hw2
m2_hi <- h2$elppd_hat + hw2

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

# 95% MC intervals
arrows(h1$S, m1_lo, h1$S, m1_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M1_A, 0.50), lwd = 2.8)
arrows(h2$S, m2_lo, h2$S, m2_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M2_A, 0.55), lwd = 2.8)

# Trajectories: M1 filled circles, M2 open circles
lines(h1$S,  h1$elppd_hat, col = COL_M1_A, lwd = 3.5)
points(h1$S, h1$elppd_hat, col = COL_M1_A, pch = 19, cex = 1.8)

lines(h2$S,  h2$elppd_hat, col = COL_M2_A, lwd = 3.5)
points(h2$S, h2$elppd_hat, col = COL_M2_A, pch = 1,  cex = 1.8, lwd = 2.2)

sep_legend_txt <- if (!is.na(sep_S)) {
  bquote("Separation (" * italic(S)[sep] == .(format(sep_S, big.mark = ",", scientific = FALSE)) * ")")
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

# Panel (b): replicated termination points (S_term, ELPPD_hat)
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

# eps = 0.05
points(t05_s_M2, t05_y_M2,
       col = adjustcolor(COL_05_M2, 0.75), pch = 19, cex = 1.05)
points(t05_s_M1, t05_y_M1,
       col = adjustcolor(COL_05_M1, 0.70), pch = 19, cex = 1.05)

# eps = 0.025
points(t025_s_M2, t025_y_M2,
       col = adjustcolor(COL_025_M2, 0.65), pch = 19, cex = 1.05)
points(t025_s_M1, t025_y_M1,
       col = adjustcolor(COL_025_M1, 0.70), pch = 19, cex = 1.05)

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

box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file1))
