## Subset selection in Bayesian linear regression (Section 6.1): non-replicated study only
## elppd_tracking_linreg.pdf: single-chain ELPPD tracking with 95% MC intervals, M1 (true) vs M2

rm(list = ls())

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

z_crit <- qnorm(0.975)

load_env <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e
}

e1 <- load_env("M1_outputs.RData")
e2 <- load_env("M2_outputs.RData")
h1 <- e1$history
h2 <- e2$history
elpd_M1 <- e1$elpd_M1   # exact LOO reference
elpd_M2 <- e1$elpd_M2

# asymp_se is stored as the MC variance in this study, hence sqrt()
hw1 <- z_crit * sqrt(pmax(h1$asymp_se, 0))
hw2 <- z_crit * sqrt(pmax(h2$asymp_se, 0))

COL_M1 <- "#08519C"   # deep royal blue
COL_M2 <- "#D73027"   # deep crimson red

S_max    <- max(h1$S, h2$S)
s_ticks  <- pretty(c(0, S_max), n = 10)
s_ticks  <- s_ticks[s_ticks > 0 & s_ticks <= S_max]
s_labels <- paste0(s_ticks / 1000, "K")

y_all   <- c(h1$elppd_hat - hw1, h1$elppd_hat + hw1,
             h2$elppd_hat - hw2, h2$elppd_hat + hw2, elpd_M1, elpd_M2)
y_rng   <- range(y_all, finite = TRUE)
y_rng   <- y_rng + c(0, 0.30) * diff(y_rng)   # headroom for the top-right legend
y_ticks <- pretty(y_rng, n = 7)

pdf_file <- "elppd_tracking_linreg.pdf"
pdf(pdf_file, width = 10, height = 8, useDingbats = FALSE)
par(mar = c(6.5, 8.5, 2.5, 2.0), oma = rep(0.5, 4), mgp = c(5.2, 1.5, 0),
    cex.lab = 2.0, cex.axis = 1.6, las = 1)

plot(NA, xlim = range(c(h1$S, h2$S)), ylim = range(y_ticks),
     xlab = expression("Iterations (" * italic(S) * ")"),
     ylab = expression(widehat(ELPPD)(italic(S))),
     xaxt = "n", yaxt = "n")
axis(1, at = s_ticks, labels = s_labels)
axis(2, at = y_ticks, labels = y_ticks)
abline(h = y_ticks, col = "grey92", lwd = 0.8)
abline(v = s_ticks, col = "grey92", lwd = 0.8)

# Exact LOO reference lines
abline(h = elpd_M1, col = COL_M1, lty = 2, lwd = 3.2)
abline(h = elpd_M2, col = COL_M2, lty = 2, lwd = 3.2)

# 95% MC intervals
arrows(h1$S, h1$elppd_hat - hw1, h1$S, h1$elppd_hat + hw1, code = 3, angle = 90,
       length = 0.06, col = adjustcolor(COL_M1, 0.50), lwd = 2.4)
arrows(h2$S, h2$elppd_hat - hw2, h2$S, h2$elppd_hat + hw2, code = 3, angle = 90,
       length = 0.06, col = adjustcolor(COL_M2, 0.55), lwd = 2.4)

# Trajectories: M1 filled circles, M2 open circles
lines(h1$S, h1$elppd_hat, col = COL_M1, lwd = 3)
points(h1$S, h1$elppd_hat, col = COL_M1, pch = 19, cex = 1.4)
lines(h2$S, h2$elppd_hat, col = COL_M2, lwd = 3)
points(h2$S, h2$elppd_hat, col = COL_M2, pch = 1, cex = 1.4, lwd = 2)

legend("topright", inset = c(0.025, 0.025),
       legend = c(expression(italic(M)[1] ~ "(True)"),
                  expression(italic(M)[2] ~ "(Misspecified)"),
                  expression(ELPPD(italic(M)[1]) ~ "(Exact)"),
                  expression(ELPPD(italic(M)[2]) ~ "(Exact)")),
       col = c(COL_M1, COL_M2, COL_M1, COL_M2), pch = c(19, 1, NA, NA),
       lty = c(1, 1, 2, 2), lwd = 3, pt.cex = 1.4, pt.lwd = c(NA, 2, NA, NA),
       bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
       cex = 1.3, y.intersp = 1.15)
box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file))
