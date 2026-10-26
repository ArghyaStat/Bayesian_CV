## Link misspecification in binary regression (Section 7.1): M1 (logistic) vs M2 (cloglog), n = 452, d = 258
## Output: elpd_mcse_arrhythmia.pdf (ELPPD tracking with 95% MC intervals and S_sep)

rm(list = ls())
invisible(gc())

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

if (!file.exists("M1_outputs.RData")) stop("M1_outputs.RData not found. Please run run_M1.R first.")
if (!file.exists("M2_outputs.RData")) stop("M2_outputs.RData not found. Please run run_M2.R first.")

load("M1_outputs.RData")
h1      <- history
d_M1    <- d
n_val   <- n
link_M1 <- if (exists("link")) link else "logit"
tau0_M1 <- if (exists("tau0")) tau0 else NA_real_
p0_M1   <- if (exists("p0")) p0 else NA_real_
rm(history)

load("M2_outputs.RData")
h2      <- history
d_M2    <- d
link_M2 <- if (exists("link")) link else "cloglog"
tau0_M2 <- if (exists("tau0")) tau0 else NA_real_
p0_M2   <- if (exists("p0")) p0 else NA_real_
rm(history)

S_vals <- h1$S
S_star <- min(S_vals)
S_max  <- max(S_vals)

## Summary
cat(sprintf("Model M1 (%s) : d = %d (D = %d shrunk), n = %d, tau0 = %.5f (p0 = %.1f)\n",
            toupper(link_M1), d_M1, d_M1 - 1, n_val, tau0_M1, p0_M1))
cat(sprintf("Model M2 (%s): d = %d (D = %d shrunk), n = %d, tau0 = %.5f (p0 = %.1f)\n",
            toupper(link_M2), d_M2, d_M2 - 1, n_val, tau0_M2, p0_M2))
cat(sprintf("Final ELPPD:\n  M1 (Logit)  : %.3f (SE = %.4f, ESS = %.1f, minESS_xi = %.1f)\n  M2 (Cloglog): %.3f (SE = %.4f, ESS = %.1f, minESS_xi = %.1f)\n  Delta (M2 - M1): %+.3f\n",
            tail(h1$elppd_hat, 1), tail(h1$asymp_se, 1),
            tail(h1$ess_elppd, 1), tail(h1$minESS_xi, 1),
            tail(h2$elppd_hat, 1), tail(h2$asymp_se, 1),
            tail(h2$ess_elppd, 1), tail(h2$minESS_xi, 1),
            tail(h2$elppd_hat, 1) - tail(h1$elppd_hat, 1)))

## Separation point S_sep
z_crit <- qnorm(0.975)

m1_lo <- h1$elppd_hat - z_crit * h1$asymp_se
m1_hi <- h1$elppd_hat + z_crit * h1$asymp_se
m2_lo <- h2$elppd_hat - z_crit * h2$asymp_se
m2_hi <- h2$elppd_hat + z_crit * h2$asymp_se

compute_S_sep <- function(S, lo1, hi1, lo2, hi2) {
  is_sep <- (lo1 > hi2) | (lo2 > hi1)
  if (length(is_sep) > 0 && is_sep[length(is_sep)]) {
    fi <- which(!is_sep)
    k_sep <- if (length(fi) == 0) 1L else max(fi) + 1L
    return(S[k_sep])
  }
  NA_real_
}

sep_S <- compute_S_sep(S_vals, m1_lo, m1_hi, m2_lo, m2_hi)
if (!is.na(sep_S)) {
  cat(sprintf("Data-adaptive separation attained at S_sep = %s iterations!\n",
              format(sep_S, big.mark = ",", scientific = FALSE)))
  k_sep <- which(S_vals == sep_S)
  cat(sprintf("  At S = %s: M1 = %.2f [%.2f, %.2f] vs M2 = %.2f [%.2f, %.2f]\n\n",
              format(sep_S, big.mark = ",", scientific = FALSE),
              h1$elppd_hat[k_sep], m1_lo[k_sep], m1_hi[k_sep],
              h2$elppd_hat[k_sep], m2_lo[k_sep], m2_hi[k_sep]))
} else {
  cat("No complete non-overlapping separation attained through S_max.\n\n")
}

## Plot styling and axes
COL_M1  <- "#08519C"
COL_M2  <- "#D73027"
COL_SEP <- "#006837"

xlim_common <- range(S_vals)
s_ticks     <- pretty(S_vals, n = 5)
s_ticks     <- s_ticks[s_ticks >= S_star & s_ticks <= S_max]
s_labels    <- if (S_max >= 1e3) paste0(s_ticks / 1000, "K") else as.character(s_ticks)

elppd_min <- min(c(h1$elppd_hat, h2$elppd_hat, m1_lo, m2_lo), na.rm = TRUE)
elppd_max <- max(c(h1$elppd_hat, h2$elppd_hat, m1_hi, m2_hi), na.rm = TRUE)
y_range   <- elppd_max - elppd_min
if (y_range < 1e-5) y_range <- 1.0

y_ticks_elppd <- pretty(c(elppd_min - 0.08 * y_range, elppd_max + 0.25 * y_range), n = 7)
ylim_elppd    <- range(y_ticks_elppd)

set_par <- function() {
  par(
    mar      = c(6.5, 8.5, 2.5, 2.0),
    oma      = c(0.5, 0.5, 0.5, 0.5),
    mgp      = c(5.2, 1.5, 0),
    cex.lab  = 2.6,
    cex.axis = 2.1,
    las      = 1
  )
}

draw_frame <- function(ylim, y_ticks, y_labels, ylab) {
  plot(NA, xlim = xlim_common, ylim = ylim,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = ylab, main = "", xaxt = "n", yaxt = "n")
  axis(1, at = s_ticks, labels = s_labels, gap.axis = 0)
  axis(2, at = y_ticks, labels = y_labels)
  abline(h = y_ticks, col = "grey92", lty = 1, lwd = 0.8)
  abline(v = s_ticks, col = "grey92", lty = 1, lwd = 0.8)
}

render_elppd_panel <- function() {
  draw_frame(ylim_elppd, y_ticks_elppd, y_ticks_elppd,
             expression(widehat(ELPPD)(italic(S))))

  if (!is.na(sep_S)) abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.2)

  arrows(h1$S, m1_lo, h1$S, pmin(m1_hi, ylim_elppd[2]), code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M1, 0.50), lwd = 2.4)
  arrows(h2$S, m2_lo, h2$S, pmin(m2_hi, ylim_elppd[2]), code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M2, 0.55), lwd = 2.4)

  lines(h1$S,  h1$elppd_hat, col = COL_M1, lwd = 3.2)
  points(h1$S, h1$elppd_hat, col = COL_M1, pch = 19, cex = 1.6)
  lines(h2$S,  h2$elppd_hat, col = COL_M2, lwd = 3.2)
  points(h2$S, h2$elppd_hat, col = COL_M2, pch = 1, cex = 1.6, lwd = 2.0)

  sep_txt <- if (!is.na(sep_S)) {
    bquote("Separation (" * italic(S)[sep] == .(format(sep_S, big.mark = ",", scientific = FALSE)) * ")")
  } else {
    expression("No separation")
  }

  legend("topright", inset = c(0.025, 0.025),
         legend = c(
           as.expression(bquote(italic(M)[1] ~ "(Logistic link)")),
           as.expression(bquote(italic(M)[2] ~ "(Cloglog link)")),
           sep_txt
         ),
         col = c(COL_M1, COL_M2, COL_SEP),
         pch = c(19, 1, NA),
         lty = c(1, 1, 4),
         lwd = c(3.0, 3.0, 3.2),
         pt.cex = c(1.6, 1.6, NA),
         pt.lwd = c(NA, 2.0, NA),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.35, y.intersp = 1.15)
  box(lwd = 1.8)
}

## Export
pdf("elpd_mcse_arrhythmia.pdf", width = 10, height = 9.5, useDingbats = FALSE)
set_par()
render_elppd_panel()
invisible(dev.off())
cat("Saved: elpd_mcse_arrhythmia.pdf\n")
