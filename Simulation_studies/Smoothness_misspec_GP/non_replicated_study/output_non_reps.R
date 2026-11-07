## Smoothness misspecification in GP regression (Section 6.2): non-replicated study only
## eljpd_tracking_smoothness.pdf: single-chain 10-fold ELJPD tracking with 95% MC intervals,
##                                M1 (nu = 0.5, true) vs M2 (nu = 1)

rm(list = ls())

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

z_crit <- qnorm(0.975)

load_env <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e
}

nr <- list(M1 = load_env("M1_outputs.RData")$history,
           M2 = load_env("M2_outputs.RData")$history)

# Reference ELJPD: exact if available, else brute force, else none
ref_label <- NULL
ref       <- NULL
if (file.exists("exact_reference.RData")) {
  ex <- load_env("exact_reference.RData")$exact
  ref <- c(M1 = ex$M1$eljpd, M2 = ex$M2$eljpd); ref_label <- "Exact"
} else if (file.exists("bruteforce_reference.RData")) {
  bf <- load_env("bruteforce_reference.RData")$bf
  ref <- c(M1 = bf$M1$eljpd, M2 = bf$M2$eljpd); ref_label <- "Brute force"
}

# asymp_se is the MC standard error in this study
ci <- lapply(nr, function(h) cbind(lo = h$eljpd_hat - z_crit * h$asymp_se,
                                   hi = h$eljpd_hat + z_crit * h$asymp_se))

COL_M1 <- "#08519C"   # deep royal blue
COL_M2 <- "#D73027"   # deep crimson red

# Equally spaced round ticks on [0, S_max]
k_ticks <- function(S_max, n_int = 5) {
  base <- 10^floor(log10(S_max / n_int))
  cand <- c(1, 2, 2.5, 3, 5, 10) * base
  cand <- cand[abs(S_max / cand - round(S_max / cand)) < 1e-8]
  step <- cand[which.min(abs(cand - S_max / n_int))]
  t <- seq(0, S_max, by = step)
  list(at = t, labels = ifelse(t == 0, "0", paste0(t / 1000, "K")))
}

xlim <- c(0, max(nr$M1$S, nr$M2$S))
ylim <- range(unlist(ci), ref)
ylim <- ylim + c(-0.05, 0.25) * diff(ylim)
xt <- k_ticks(xlim[2])
yt <- pretty(ylim)

pdf_file <- "eljpd_tracking_smoothness.pdf"
pdf(pdf_file, width = 10, height = 8, useDingbats = FALSE)
par(mar = c(6.5, 10, 2.5, 2.0), oma = rep(0.5, 4), mgp = c(5.2, 1.5, 0),
    cex.lab = 2.0, cex.axis = 1.6, las = 1)

plot(NA, xlim = xlim, ylim = ylim, xaxt = "n",
     xlab = expression("Iterations (" * italic(S) * ")"), ylab = "")
title(ylab = expression(widehat(ELJPD)(italic(S))), line = 7)
axis(1, at = xt$at, labels = xt$labels)
abline(h = yt, col = "grey92", lwd = 0.8)
abline(v = xt$at, col = "grey92", lwd = 0.8)

if (!is.null(ref)) {
  abline(h = ref["M1"], col = COL_M1, lty = 2, lwd = 3.2)
  abline(h = ref["M2"], col = COL_M2, lty = 2, lwd = 3.2)
}

arrows(nr$M1$S, ci$M1[, "lo"], nr$M1$S, ci$M1[, "hi"], code = 3, angle = 90,
       length = 0.06, col = adjustcolor(COL_M1, 0.50), lwd = 2.4)
arrows(nr$M2$S, ci$M2[, "lo"], nr$M2$S, ci$M2[, "hi"], code = 3, angle = 90,
       length = 0.06, col = adjustcolor(COL_M2, 0.55), lwd = 2.4)
lines(nr$M1$S, nr$M1$eljpd_hat, col = COL_M1, lwd = 3)
points(nr$M1$S, nr$M1$eljpd_hat, col = COL_M1, pch = 19, cex = 1.4)
lines(nr$M2$S, nr$M2$eljpd_hat, col = COL_M2, lwd = 3)
points(nr$M2$S, nr$M2$eljpd_hat, col = COL_M2, pch = 1, cex = 1.4, lwd = 2)

leg <- list(legend = c(expression(italic(M)[1] ~ (nu == 0.5) ~ "(True)"),
                       expression(italic(M)[2] ~ (nu == 1) ~ "(Misspecified)")),
            col = c(COL_M1, COL_M2), pch = c(19, 1), lty = c(1, 1),
            lwd = c(3, 3), pt.lwd = c(NA, 2))
if (!is.null(ref)) {
  leg$legend <- c(leg$legend,
                  as.expression(bquote(ELJPD(italic(M)[1]) ~ "(" * .(ref_label) * ")")),
                  as.expression(bquote(ELJPD(italic(M)[2]) ~ "(" * .(ref_label) * ")")))
  leg$col <- c(leg$col, COL_M1, COL_M2); leg$pch <- c(leg$pch, NA, NA)
  leg$lty <- c(leg$lty, 2, 2); leg$lwd <- c(leg$lwd, 3.2, 3.2)
  leg$pt.lwd <- c(leg$pt.lwd, NA, NA)
}
do.call(legend, c(list("topright", inset = c(0.025, 0.025), pt.cex = 1.4, bty = "o",
                       bg = adjustcolor("white", 0.94), box.col = "grey75",
                       cex = 1.3, y.intersp = 1.15), leg))
box(lwd = 1.8)

dev.off()
cat(sprintf("Saved %s\n", pdf_file))
