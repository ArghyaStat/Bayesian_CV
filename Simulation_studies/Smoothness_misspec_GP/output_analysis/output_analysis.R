## Smoothness misspecification in GP regression (Section 6.2): M1 (nu = 0.5, true) vs M2 (nu = 1), 10-fold ELJPD
## elpd_mcse_smoothness.pdf: (a) single-chain ELJPD tracking with 95% MC intervals and S_sep;
##                           (b) replicated termination points (S_term, ELJPD_hat) for eps = 0.05, 0.025

rm(list = ls())

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

nonrep_dir <- "../non_replicated_study"
ref_dir    <- nonrep_dir
rep_dir    <- "../replicated_study"

eps_vals   <- c(0.05, 0.025)
ess_target <- c(6147, 24585)
z_crit     <- qnorm(0.975)

load_env <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e
}

nr <- list(M1 = load_env(file.path(nonrep_dir, "M1_outputs.RData"))$history,
           M2 = load_env(file.path(nonrep_dir, "M2_outputs.RData"))$history)
rp <- list(M1 = load_env(file.path(rep_dir, "M1_outputs_reps.RData"))$results_reps,
           M2 = load_env(file.path(rep_dir, "M2_outputs_reps.RData"))$results_reps)

exact_file <- file.path(ref_dir, "exact_reference.RData")
bf_file    <- file.path(ref_dir, "bruteforce_reference.RData")
# Reference ELJPD: exact if available, else brute force
ref_label  <- if (file.exists(exact_file)) "Exact" else "Brute force"
ref <- if (file.exists(exact_file)) {
  ex <- load_env(exact_file)$exact
  c(M1 = ex$M1$eljpd, M2 = ex$M2$eljpd)
} else if (file.exists(bf_file)) {
  bf <- load_env(bf_file)$bf
  c(M1 = bf$M1$eljpd, M2 = bf$M2$eljpd)
} else NULL

# First checkpoint with ESS_ELJPD >= target (S_max if never reached)
get_term <- function(h, target) {
  k <- which(h$ess_elpd >= target)[1]
  reached <- !is.na(k)
  if (!reached) k <- length(h$S)
  data.frame(S_term = h$S[k], ELJPD_hat = h$eljpd_hat[k],
             asymp_se = h$asymp_se[k], reached = reached)
}

# Earliest checkpoint after which the 95% intervals of M1 and M2 stay disjoint
compute_S_sep <- function(h1, h2) {
  lo1 <- h1$eljpd_hat - z_crit * h1$asymp_se
  hi1 <- h1$eljpd_hat + z_crit * h1$asymp_se
  lo2 <- h2$eljpd_hat - z_crit * h2$asymp_se
  hi2 <- h2$eljpd_hat + z_crit * h2$asymp_se
  is_sep <- lo1 > hi2 | lo2 > hi1
  if (!is_sep[length(is_sep)]) return(NA_real_)
  fi <- which(!is_sep)
  h1$S[if (length(fi) == 0) 1L else max(fi) + 1L]
}

term_nr <- lapply(nr, function(h) lapply(ess_target, get_term, h = h))
term_rp <- lapply(rp, function(res) lapply(ess_target, function(tg)
  do.call(rbind, lapply(res, function(x) get_term(x$history, tg)))))
S_sep <- compute_S_sep(nr$M1, nr$M2)

# Summary of replicated termination points
summ <- do.call(rbind, lapply(names(term_rp), function(mod) {
  do.call(rbind, lapply(seq_along(eps_vals), function(j) {
    d <- term_rp[[mod]][[j]]
    out <- data.frame(model = mod, eps = eps_vals[j], ESS_target = ess_target[j],
                      median_S_term = median(d$S_term),
                      mean_ELJPD = mean(d$ELJPD_hat), sd_ELJPD = sd(d$ELJPD_hat),
                      mean_mcse = mean(d$asymp_se), not_reached = mean(!d$reached))
    if (!is.null(ref)) {
      lo <- d$ELJPD_hat - z_crit * d$asymp_se
      hi <- d$ELJPD_hat + z_crit * d$asymp_se
      out$coverage <- mean(lo <= ref[mod] & ref[mod] <= hi)
    }
    out
  }))
}))
rank_ok <- sapply(seq_along(eps_vals), function(j)
  mean(term_rp$M1[[j]]$ELJPD_hat > term_rp$M2[[j]]$ELJPD_hat))

if (!is.null(ref))
  cat(sprintf("%s ELJPD  M1: %.4f  M2: %.4f  difference: %.4f\n", ref_label,
              ref["M1"], ref["M2"], ref["M1"] - ref["M2"]))
cat(sprintf("Non-replicated separation point S_sep: %s\n", format(S_sep, big.mark = ",", scientific = FALSE)))
print(summ, digits = 4, row.names = FALSE)
cat(sprintf("P(ELJPD_hat M1 > M2), eps = %s: %.3f\n", eps_vals, rank_ok), sep = "")

COL_M1_A <- "#08519C"   # deep royal blue
COL_M2_A <- "#D73027"   # deep crimson red
COL_SEP  <- "#006837"   # forest green
COL_05   <- c(M1 = "#D73027", M2 = "#D4AC0D")   # crimson / golden-yellow
COL_025  <- c(M1 = "#08519C", M2 = "#6BAED6")   # royal blue / sky blue

# Equally spaced round ticks on [0, S_max], e.g. 30K steps for S_max = 150K
k_ticks <- function(S_max, n_int = 5) {
  base  <- 10^floor(log10(S_max / n_int))
  cand  <- c(1, 2, 2.5, 3, 5, 10) * base
  cand  <- cand[abs(S_max / cand - round(S_max / cand)) < 1e-8]
  step  <- cand[which.min(abs(cand - S_max / n_int))]
  t <- seq(0, S_max, by = step)
  list(at = t, labels = ifelse(t == 0, "0", paste0(t / 1000, "K")))
}

ci <- lapply(nr, function(h) cbind(lo = h$eljpd_hat - z_crit * h$asymp_se,
                                   hi = h$eljpd_hat + z_crit * h$asymp_se))

pdf("elpd_mcse_smoothness.pdf", width = 20, height = 9.5, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(6.5, 11, 2.5, 2.0), oma = rep(0.5, 4),
    mgp = c(5.2, 1.5, 0), cex.lab = 2.6, cex.axis = 2.1, las = 1)

# Panel (a): single-chain ELJPD tracking and separation point
xlim1 <- c(0, max(nr$M1$S))
ylim1 <- range(unlist(ci), ref)
ylim1 <- ylim1 + c(-0.05, 0.25) * diff(ylim1)
xt1 <- k_ticks(xlim1[2])
yt1 <- pretty(ylim1)

plot(NA, xlim = xlim1, ylim = ylim1, xaxt = "n",
     xlab = expression("Iterations (" * italic(S) * ")"), ylab = "")
title(ylab = expression(widehat(ELJPD)(italic(S))), line = 8)
axis(1, at = xt1$at, labels = xt1$labels)
abline(h = yt1, col = "grey92", lwd = 0.8)
abline(v = xt1$at, col = "grey92", lwd = 0.8)

if (!is.null(ref)) {
  abline(h = ref["M1"], col = COL_M1_A, lty = 2, lwd = 3.2)
  abline(h = ref["M2"], col = COL_M2_A, lty = 2, lwd = 3.2)
}
if (!is.na(S_sep)) abline(v = S_sep, col = COL_SEP, lty = 4, lwd = 3.5)

arrows(nr$M1$S, ci$M1[, "lo"], nr$M1$S, ci$M1[, "hi"], code = 3, angle = 90,
       length = 0.08, col = adjustcolor(COL_M1_A, 0.50), lwd = 2.8)
arrows(nr$M2$S, ci$M2[, "lo"], nr$M2$S, ci$M2[, "hi"], code = 3, angle = 90,
       length = 0.08, col = adjustcolor(COL_M2_A, 0.55), lwd = 2.8)
lines(nr$M1$S, nr$M1$eljpd_hat, col = COL_M1_A, lwd = 3.5)
points(nr$M1$S, nr$M1$eljpd_hat, col = COL_M1_A, pch = 19, cex = 1.8)
lines(nr$M2$S, nr$M2$eljpd_hat, col = COL_M2_A, lwd = 3.5)
points(nr$M2$S, nr$M2$eljpd_hat, col = COL_M2_A, pch = 1, cex = 1.8, lwd = 2.2)

sep_txt <- if (!is.na(S_sep)) {
  bquote("Separation (" * italic(S)[sep] == .(format(S_sep, big.mark = ",", scientific = FALSE)) * ")")
} else expression("No separation")

leg_a <- list(legend = c(expression(italic(M)[1] ~ (nu == 0.5) ~ "(True)"),
                         expression(italic(M)[2] ~ (nu == 1) ~ "(Misspecified)")),
              col = c(COL_M1_A, COL_M2_A), pch = c(19, 1), lty = c(1, 1),
              lwd = c(3.2, 3.2), pt.cex = c(1.8, 1.8), pt.lwd = c(NA, 2.2))
if (!is.null(ref)) {
  leg_a$legend <- c(leg_a$legend,
                    as.expression(bquote(ELJPD(italic(M)[1]) ~ "(" * .(ref_label) * ")")),
                    as.expression(bquote(ELJPD(italic(M)[2]) ~ "(" * .(ref_label) * ")")))
  leg_a$col <- c(leg_a$col, COL_M1_A, COL_M2_A); leg_a$pch <- c(leg_a$pch, NA, NA)
  leg_a$lty <- c(leg_a$lty, 2, 2); leg_a$lwd <- c(leg_a$lwd, 3.2, 3.2)
  leg_a$pt.cex <- c(leg_a$pt.cex, NA, NA); leg_a$pt.lwd <- c(leg_a$pt.lwd, NA, NA)
}
leg_a$legend <- c(leg_a$legend, sep_txt)
leg_a$col <- c(leg_a$col, COL_SEP); leg_a$pch <- c(leg_a$pch, NA)
leg_a$lty <- c(leg_a$lty, 4); leg_a$lwd <- c(leg_a$lwd, 3.5)
leg_a$pt.cex <- c(leg_a$pt.cex, NA); leg_a$pt.lwd <- c(leg_a$pt.lwd, NA)
do.call(legend, c(list("topright", inset = c(0.025, 0.025), bty = "o",
                       bg = adjustcolor("white", 0.94), box.col = "grey75",
                       cex = 1.55, y.intersp = 1.15), leg_a))
box(lwd = 1.8)

# Panel (b): replicated termination distributions
s_all <- unlist(lapply(term_rp, function(x) lapply(x, `[[`, "S_term")))
y_all <- unlist(lapply(term_rp, function(x) lapply(x, `[[`, "ELJPD_hat")))
xlim2 <- c(0, max(rp$M1[[1]]$history$S))
ylim2 <- range(y_all, ref)
ylim2 <- ylim2 + c(-0.05, 0.30) * diff(ylim2)
xt2 <- k_ticks(xlim2[2])
yt2 <- pretty(ylim2)
jit <- 0.15 * diff(rp$M1[[1]]$history$S[1:2])

plot(NA, xlim = xlim2, ylim = ylim2, xaxt = "n",
     xlab = expression("Iterations for termination (" * italic(S)[term] * ")"), ylab = "")
title(ylab = expression(widehat(ELJPD)(italic(S)[term])), line = 8)
axis(1, at = xt2$at, labels = xt2$labels)
abline(h = yt2, col = "grey92", lwd = 0.8)
abline(v = xt2$at, col = "grey92", lwd = 0.8)
if (!is.null(ref)) abline(h = ref, col = "grey35", lty = 2, lwd = 3.0)

points(jitter(term_rp$M2[[1]]$S_term, amount = jit), term_rp$M2[[1]]$ELJPD_hat,
       col = adjustcolor(COL_05["M2"], 0.75), pch = 19, cex = 1.05)
points(jitter(term_rp$M1[[1]]$S_term, amount = jit), term_rp$M1[[1]]$ELJPD_hat,
       col = adjustcolor(COL_05["M1"], 0.70), pch = 19, cex = 1.05)
points(jitter(term_rp$M2[[2]]$S_term, amount = jit), term_rp$M2[[2]]$ELJPD_hat,
       col = adjustcolor(COL_025["M2"], 0.65), pch = 19, cex = 1.05)
points(jitter(term_rp$M1[[2]]$S_term, amount = jit), term_rp$M1[[2]]$ELJPD_hat,
       col = adjustcolor(COL_025["M1"], 0.70), pch = 19, cex = 1.05)

leg_b <- c(expression(epsilon == 0.05 * ":" ~ italic(M)[1] ~ "(True)"),
           expression(epsilon == 0.05 * ":" ~ italic(M)[2] ~ "(Misspecified)"),
           expression(epsilon == 0.025 * ":" ~ italic(M)[1] ~ "(True)"),
           expression(epsilon == 0.025 * ":" ~ italic(M)[2] ~ "(Misspecified)"))
col_b <- c(COL_05["M1"], COL_05["M2"], COL_025["M1"], COL_025["M2"])
if (!is.null(ref)) {
  legend("topright", inset = c(0.025, 0.025),
         legend = c(leg_b, as.expression(bquote(.(ref_label) ~ ELJPD ~ "(Reference)"))),
         col = c(col_b, "grey35"), pch = c(19, 19, 19, 19, NA), lty = c(NA, NA, NA, NA, 2),
         lwd = c(NA, NA, NA, NA, 3.0), pt.cex = c(1.5, 1.5, 1.5, 1.5, NA),
         bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.35, y.intersp = 1.15)
} else {
  legend("topright", inset = c(0.025, 0.025), legend = leg_b, col = col_b, pch = 19,
         pt.cex = 1.5, bty = "o", bg = adjustcolor("white", 0.94), box.col = "grey75",
         cex = 1.35, y.intersp = 1.15)
}
box(lwd = 1.8)

dev.off()
cat("Saved elpd_mcse_smoothness.pdf\n")
