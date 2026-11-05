## Subset selection in Bayesian linear regression (Section 6.1): M1 (true) vs M2 (one predictor dropped)
## elpd_mcse_linreg.pdf: (a) single-chain ELPPD tracking with 95% MC intervals;
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

# Single-chain runs
load(file.path(nonrep_dir, "M1_outputs.RData"))
nr_M1_history <- history
nr_elpd_M1    <- elpd_M1
nr_elpd_M2    <- elpd_M2
rm(history, acc_mix, elpd_diff_val)

load(file.path(nonrep_dir, "M2_outputs.RData"))
nr_M2_history <- history
rm(history, acc_mix, elpd_diff_val, R2_dropped, tau2_dropped, dropped_col)

# Replicated runs
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
  ess    <- if (!is.null(h$ess_elppd)) h$ess_elppd else h$ess
  idx_ok <- which(ess >= ess_target)
  k      <- if (length(idx_ok) > 0) idx_ok[1] else length(h$S)
  list(S_term    = h$S[k],
       ELPPD_hat = h$elppd_hat[k],
       asymp_se  = h$asymp_se[k],
       k         = k)
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
COL_05_M1 <- "#D73027"   # deep crimson red
COL_05_M2 <- "#D4AC0D"   # rich golden-yellow
COL_025_M1 <- "#08519C"   # deep royal blue
COL_025_M2 <- "#6BAED6"   # soft sky-blue tint

# Axis ranges derived from the data (shared y-range across both panels)
S_all <- c(nr_M1_history$S, nr_M2_history$S,
           unlist(lapply(c(rep_M1, rep_M2), function(x) x$history$S)))
s_ticks_uniform  <- pretty(c(0, max(S_all)), n = 10)
s_ticks_uniform  <- s_ticks_uniform[s_ticks_uniform > 0 & s_ticks_uniform <= max(S_all)]
s_labels_uniform <- paste0(s_ticks_uniform / 1000, "K")

y_all <- c(nr_M1_history$elppd_hat - z_crit * sqrt(pmax(nr_M1_history$asymp_se, 0)),
           nr_M1_history$elppd_hat + z_crit * sqrt(pmax(nr_M1_history$asymp_se, 0)),
           nr_M2_history$elppd_hat - z_crit * sqrt(pmax(nr_M2_history$asymp_se, 0)),
           nr_M2_history$elppd_hat + z_crit * sqrt(pmax(nr_M2_history$asymp_se, 0)),
           t05_y_M1, t05_y_M2, t025_y_M1, t025_y_M2,
           nr_elpd_M1, nr_elpd_M2, rep_elpd_M1, rep_elpd_M2)
y_ticks_all <- pretty(range(y_all, finite = TRUE), n = 7)
ylim_all    <- range(y_ticks_all)

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

# Panel (a): single-chain ELPPD tracking
h1 <- nr_M1_history
h2 <- nr_M2_history

hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))

m1_lo <- h1$elppd_hat - hw1
m1_hi <- h1$elppd_hat + hw1
m2_lo <- h2$elppd_hat - hw2
m2_hi <- h2$elppd_hat + hw2

ylim1 <- ylim_all
xlim1 <- range(h1$S)

plot(NA, xlim = xlim1, ylim = ylim1,
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ELPPD)(italic(S))),
     main = "", xaxt = "n", yaxt = "n")

s_ticks1  <- s_ticks_uniform
s_labels1 <- s_labels_uniform
axis(1, at = s_ticks1, labels = s_labels1)

y_ticks1 <- y_ticks_all
axis(2, at = y_ticks1, labels = y_ticks1)

# Subtle background grid
abline(h = y_ticks1, col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks1, col = "grey92", lty = 1, lwd = 0.8)

# Exact LOO reference lines
abline(h = nr_elpd_M1, col = COL_M1_A, lty = 2, lwd = 3.2)
abline(h = nr_elpd_M2, col = COL_M2_A, lty = 2, lwd = 3.2)

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

legend("topright", inset = c(0.025, 0.025),
       legend = c(
         expression(italic(M)[1] ~ "(True)"),
         expression(italic(M)[2] ~ "(Misspecified)"),
         expression(ELPPD(italic(M)[1]) ~ "(Exact)"),
         expression(ELPPD(italic(M)[2]) ~ "(Exact)")
       ),
       col = c(COL_M1_A, COL_M2_A, COL_M1_A, COL_M2_A),
       pch = c(19, 1, NA, NA),
       lty = c(1, 1, 2, 2),
       lwd = c(3.2, 3.2, 3.2, 3.2),
       pt.cex = c(1.8, 1.8, NA, NA),
       pt.lwd = c(NA, 2.2, NA, NA),
       bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex = 1.55, y.intersp = 1.15)

box(lwd = 1.8)

# Panel (b): replicated termination points (S_term, ELPPD_hat)
xlim2 <- range(c(t05_s_M1, t05_s_M2, t025_s_M1, t025_s_M2,
                 rep_M1[[1]]$history$S))
ylim2 <- ylim_all

plot(NA, xlim = xlim2, ylim = ylim2,
     xlab = expression("Iterations for termination (" * italic(S)[term] * ")"),
     ylab = expression(widehat(ELPPD)(italic(S)[term])),
     main = "", xaxt = "n", yaxt = "n")

s_ticks2  <- s_ticks_uniform
s_labels2 <- s_labels_uniform
axis(1, at = s_ticks2, labels = s_labels2)

y_ticks2 <- y_ticks_all
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
