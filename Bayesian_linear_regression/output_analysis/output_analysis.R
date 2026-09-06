# High-dimensional case study: n = 500, d = 250 (Mixture IS)
# Panel (a): Single-chain ELPPD tracking with 95% CI and separation point
# Panel (b): Replicated termination distributions (S_term, ELPPD_hat)

rm(list = ls())
gc()

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

# Load non-replicated results
load("d250_n500_true_results.RData")
nr_true    <- results
nr_elpd_M1 <- elpd_M1
nr_elpd_M2 <- elpd_M2

load("d250_n500_misspec_results.RData")
nr_mis <- results

# Load replicated results (100 reps)
load("d250_n500_true_results_reps100.RData")
rep_true    <- results_reps
rep_elpd_M1 <- elpd_M1
rep_elpd_M2 <- elpd_M2

load("d250_n500_misspec_results_reps100.RData")
rep_mis <- results_reps

pdf_file <- "elppd_side_by_side_d250_n500.pdf"
pdf(pdf_file, width = 20, height = 9.5, useDingbats = FALSE)

par(
  mfrow    = c(1, 2),
  mar      = c(6.5, 7.8, 2.5, 2.0),
  oma      = c(0.5, 0.5, 0.5, 0.5),
  mgp      = c(4.8, 1.5, 0),
  cex.lab  = 2.6,
  cex.axis = 2.1,
  las      = 1
)

# Panel (a): Dynamic ELPPD tracking and point of separation

h1 <- nr_true$history
h2 <- nr_mis$history

CI_Z <- qnorm(0.975)
hw1  <- CI_Z * sqrt(pmax(h1$asymp_se, 0))
hw2  <- CI_Z * sqrt(pmax(h2$asymp_se, 0))

m1_lo <- h1$elppd_hat - hw1
m1_hi <- h1$elppd_hat + hw1
m2_lo <- h2$elppd_hat - hw2
m2_hi <- h2$elppd_hat + hw2

# Separation: earliest S after which m1_lo > m2_hi through S_max
is_sep <- (m1_lo > m2_hi)
if (is_sep[length(is_sep)]) {
  false_idx <- which(!is_sep)
  k_sep     <- if (length(false_idx) == 0) 1L else max(false_idx) + 1L
  sep_S     <- h1$S[k_sep]
} else {
  sep_S <- NA_real_
}

ylim1 <- c(-874.5, -837.0)
xlim1 <- range(h1$S)

COL_M1  <- "#08519C"
COL_M2  <- "#D73027"
COL_SEP <- "#006837"

plot(NA, xlim = xlim1, ylim = ylim1,
     xlab = expression("MCMC Draws" ~ italic(S)),
     ylab = expression(widehat(ELPPD)(italic(S))),
     main = "", xaxt = "n")

s_ticks  <- c(5e3, 25e3, 50e3, 75e3, 1e5)
s_labels <- c("5K", "25K", "50K", "75K", "100K")
axis(1, at = s_ticks, labels = s_labels)

abline(h = pretty(ylim1, n = 6), col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks,               col = "grey92", lty = 1, lwd = 0.8)

abline(h = nr_elpd_M1, col = COL_M1, lty = 2, lwd = 3.2)
abline(h = nr_elpd_M2, col = COL_M2, lty = 2, lwd = 3.2)

if (!is.na(sep_S)) {
  abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.5)
}

arrows(h1$S, m1_lo, h1$S, m1_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M1, 0.55), lwd = 2.8)
arrows(h2$S, m2_lo, h2$S, m2_hi, code = 3, angle = 90, length = 0.08,
       col = adjustcolor(COL_M2, 0.55), lwd = 2.8)

lines(h1$S,  h1$elppd_hat, col = COL_M1, lwd = 3.5)
points(h1$S, h1$elppd_hat, col = COL_M1, pch = 19, cex = 1.9)

lines(h2$S,  h2$elppd_hat, col = COL_M2, lwd = 3.5)
points(h2$S, h2$elppd_hat, col = COL_M2, pch = 17, cex = 1.9)

text(x = xlim1[1] + 0.03 * diff(xlim1), y = ylim1[2] - 0.04 * diff(ylim1),
     labels = "(a)", cex = 2.6, font = 2, adj = c(0, 1))

legend("topright", inset = c(0.025, 0.025),
       legend = c(
         expression(italic(M)[1] ~ "(True)"),
         expression(italic(M)[2] ~ "(Misspecified)"),
         expression(ELPPD(italic(M)[1]) ~ "(Exact)"),
         expression(ELPPD(italic(M)[2]) ~ "(Exact)"),
         bquote("Separation (" * italic(S)[sep] == 
                  .(format(sep_S, big.mark = ",")) * ")")
       ),
       col = c(COL_M1, COL_M2, COL_M1, COL_M2, COL_SEP),
       pch = c(19, 17, NA, NA, NA),
       lty = c(1, 1, 2, 2, 4),
       lwd = c(3.2, 3.2, 3.2, 3.2, 3.5),
       pt.cex = c(1.8, 1.8, NA, NA, NA),
       bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex = 1.65, y.intersp = 1.15)

box(lwd = 1.8)

# Panel (b): Replicated termination distributions (S_term, ELPPD_hat)

t05_s <- sapply(rep_true, function(x) x$summary_05$S_term)
t05_y <- sapply(rep_true, function(x) x$summary_05$ELPPD_hat)
t03_s <- sapply(rep_true, function(x) x$summary_03$S_term)
t03_y <- sapply(rep_true, function(x) x$summary_03$ELPPD_hat)

m05_s <- sapply(rep_mis, function(x) x$summary_05$S_term)
m05_y <- sapply(rep_mis, function(x) x$summary_05$ELPPD_hat)
m03_s <- sapply(rep_mis, function(x) x$summary_03$S_term)
m03_y <- sapply(rep_mis, function(x) x$summary_03$ELPPD_hat)

# eps = 0.05: amber; eps = 0.03: blue
COL_05_M1 <- "#E65100"
COL_05_M2 <- "#FDAE61"
COL_03_M1 <- "#08519C"
COL_03_M2 <- "#9ECAE1"

xlim2 <- c(24e3, 156e3)
ylim2 <- c(-856.2, -841.2)

plot(NA, xlim = xlim2, ylim = ylim2,
     xlab = expression("Termination Effort" ~ italic(S)[term]),
     ylab = expression(widehat(ELPPD)(italic(S)[term])),
     main = "", xaxt = "n")

s_ticks2  <- c(30e3, 60e3, 90e3, 120e3, 150e3)
s_labels2 <- c("30K", "60K", "90K", "120K", "150K")
axis(1, at = s_ticks2, labels = s_labels2)

abline(h = pretty(c(-852, -842), n = 6), col = "grey92", lty = 1, lwd = 0.8)
abline(v = s_ticks2,                     col = "grey92", lty = 1, lwd = 0.8)

abline(h = rep_elpd_M1, col = "grey35", lty = 2, lwd = 3.0)
abline(h = rep_elpd_M2, col = "grey35", lty = 2, lwd = 3.0)

set.seed(42)
jit_w <- 1800

points(jitter(t05_s, amount = jit_w), t05_y,
       col = adjustcolor(COL_05_M1, 0.88), pch = 19, cex = 1.7)
points(jitter(m05_s, amount = jit_w), m05_y,
       col = adjustcolor("#B33D00", 0.95), bg = adjustcolor(COL_05_M2, 0.75),
       pch = 24, cex = 1.7, lwd = 1.8)

points(jitter(t03_s, amount = jit_w), t03_y,
       col = adjustcolor(COL_03_M1, 0.88), pch = 19, cex = 1.7)
points(jitter(m03_s, amount = jit_w), m03_y,
       col = adjustcolor("#084081", 0.95), bg = adjustcolor(COL_03_M2, 0.75),
       pch = 24, cex = 1.7, lwd = 1.8)

text(x = xlim2[1] + 0.03 * diff(xlim2), y = ylim2[2] - 0.04 * diff(ylim2),
     labels = "(b)", cex = 2.6, font = 2, adj = c(0, 1))

legend("bottomleft", inset = c(0.025, 0.025),
       legend = c(
         expression(epsilon == 0.05 * ":" ~ italic(M)[1] ~ "(True)"),
         expression(epsilon == 0.05 * ":" ~ italic(M)[2] ~ "(Misspecified)"),
         expression(epsilon == 0.03 * ":" ~ italic(M)[1] ~ "(True)"),
         expression(epsilon == 0.03 * ":" ~ italic(M)[2] ~ "(Misspecified)"),
         expression("Exact" ~ ELPPD ~ "(Reference)")
       ),
       col = c(COL_05_M1, "#B33D00", COL_03_M1, "#084081", "grey35"),
       pt.bg = c(COL_05_M1, adjustcolor(COL_05_M2, 0.75),
                 COL_03_M1, adjustcolor(COL_03_M2, 0.75), NA),
       pch = c(19, 24, 19, 24, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, 1.8, NA, 1.8, 3.0),
       pt.cex = c(1.7, 1.6, 1.7, 1.6, NA),
       bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex = 1.55, y.intersp = 1.15)

box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file))
