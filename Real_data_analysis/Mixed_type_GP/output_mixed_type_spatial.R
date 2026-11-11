## Latent covariance misspecification in mixed-type spatial data (Section 7.2): US wildfires, N = 3,503 grid cells
## M1: joint (separable cross-covariance) GP, M2: independent GPs; output elpd_mcse_spatial.pdf

rm(list = ls())
invisible(gc())

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

## Which run to summarise:
##   ""        -> run_M1.R / run_M2.R               (M*_outputs.RData; default)
##   "_burnin" -> run_M1_burnin.R / run_M2_burnin.R (M*_outputs_burnin.RData)
## Can also be passed on the command line: Rscript output_mixed_type_spatial.R _burnin
run_tag <- ""
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1) run_tag <- cmd_args[1]

## Checkpoints of the saved history to keep, c(first, last) in sampler iterations.
## The kept checkpoints are re-indexed from 0, e.g. c(6e4, 9e4) -> 0, 1K, ..., 30K.
## Set to NULL to keep every checkpoint on its original index.
S_window <- if (run_tag == "") c(6e4, 9e4) else NULL

## Load model outputs
cat("==============================================================================\n")
cat("Real Data Analysis: Model Comparison & Convergence Diagnostics\n")
cat("==============================================================================\n\n")

out_m1 <- paste0("M1_outputs", run_tag, ".RData")
out_m2 <- paste0("M2_outputs", run_tag, ".RData")
file_m1_candidates <- c(out_m1, file.path("M1_estimation", out_m1))
file_m2_candidates <- c(out_m2, file.path("M2_estimation", out_m2))

file_m1 <- file_m1_candidates[file.exists(file_m1_candidates)][1]
file_m2 <- file_m2_candidates[file.exists(file_m2_candidates)][1]

if (is.na(file_m1) || is.na(file_m2)) {
  stop(sprintf("Output files missing! Please ensure %s and %s are present.", out_m1, out_m2))
}

env_m1 <- new.env()
env_m2 <- new.env()
load(file_m1, envir = env_m1)
load(file_m2, envir = env_m2)

h1 <- as.data.frame(env_m1$history)
h2 <- as.data.frame(env_m2$history)

extract_field <- function(df, candidates, default = NULL) {
  for (cand in candidates) {
    if (cand %in% names(df)) return(df[[cand]])
  }
  default
}

get_var <- function(env, names, default = NA_real_) {
  for (nm in names) {
    if (exists(nm, envir = env)) return(get(nm, envir = env))
  }
  default
}

## Subset the checkpoints to S_window; S_iter keeps the sampler iteration, S the new index
h1 <- h1[, setdiff(names(h1), grep("^fold_log_pred", names(h1), value = TRUE)), drop = FALSE]
h2 <- h2[, setdiff(names(h2), grep("^fold_log_pred", names(h2), value = TRUE)), drop = FALSE]
h1$S_iter <- h1$S
h2$S_iter <- h2$S
if (!is.null(S_window)) {
  h1 <- h1[h1$S >= S_window[1] & h1$S <= S_window[2], , drop = FALSE]
  h2 <- h2[h2$S >= S_window[1] & h2$S <= S_window[2], , drop = FALSE]
  if (nrow(h1) == 0 || nrow(h2) == 0) stop("No checkpoints fall inside S_window.")
  h1$S <- h1$S - S_window[1]
  h2$S <- h2$S - S_window[1]
}

h1_elpd     <- extract_field(h1, c("eljpd_hat", "elpd_hat"))
h1_asymp_se <- extract_field(h1, c("asymp_se", "se"))
h1_ess_pred <- extract_field(h1, c("ess_elpd", "ess"))
h1_min_ess  <- extract_field(h1, "minESS_xi", default = rep(NA_real_, nrow(h1)))

h2_elpd     <- extract_field(h2, c("eljpd_hat", "elpd_hat"))
h2_asymp_se <- extract_field(h2, c("asymp_se", "se"))
h2_ess_pred <- extract_field(h2, c("ess_elpd", "ess"))
h2_min_ess  <- extract_field(h2, "minESS_xi", default = rep(NA_real_, nrow(h2)))

## Each model keeps its own checkpoint grid; S is the plotted index (see S_window)
S1_vals       <- h1$S
S2_vals       <- h2$S
S_vals        <- S1_vals
S_star        <- min(S1_vals, S2_vals)
S_max         <- max(S1_vals, S2_vals)
S_burn_m1     <- get_var(env_m1, "S_burn", 0)
S_burn_m2     <- get_var(env_m2, "S_burn", 0)
from_zero     <- !is.null(S_window) || (S_burn_m1 > 0) || (S_burn_m2 > 0)
S_step        <- if (length(S_vals) > 1) diff(S_vals)[1] else NA_real_
n_checkpoints <- length(S_vals)

k_folds   <- get_var(env_m1, "k.folds", 10)
m_vecchia <- get_var(env_m1, "m", 20)

acc_m1 <- get_var(env_m1, c("acc.phi", "acc_phi"))
acc_m2 <- get_var(env_m2, c("acc.phi", "acc_phi"))

cat(sprintf("Loaded M1 from '%s' | M2 from '%s'\n", file_m1, file_m2))
cat(sprintf("Parameters: K-folds = %d, Vecchia conditioning m = %d\n", k_folds, m_vecchia))
cat(sprintf("Iterations: S_star = %s, S_max = %s, S_step = %s (%d checkpoints)\n",
            format(as.integer(S_star), big.mark = ","),
            format(as.integer(S_max), big.mark = ","),
            format(as.integer(S_step), big.mark = ","),
            n_checkpoints))
if (!is.null(S_window)) {
  cat(sprintf("Checkpoints shown: sampler iterations %s to %s, re-indexed from 0\n",
              format(as.integer(S_window[1]), big.mark = ","),
              format(as.integer(S_window[2]), big.mark = ",")))
}
cat(sprintf("Metropolis acceptance rate (phi): M1 = %.4f, M2 = %.4f\n\n", acc_m1, acc_m2))

## Summary at S_max
last_k1 <- length(S1_vals)
last_k2 <- length(S2_vals)
S_last  <- c(h1$S_iter[last_k1], h2$S_iter[last_k2])   # sampler iteration of the last checkpoint shown

summary_df <- data.frame(
  Model           = c("M1 (Joint)", "M2 (Separate)"),
  Specification   = c("Separable Cross-Covariance GP", "Independent Univariate GPs"),
  S_max           = format(as.integer(S_last), big.mark = ","),
  ELPD_hat        = round(c(h1_elpd[last_k1], h2_elpd[last_k2]), 2),
  MCSE            = round(c(h1_asymp_se[last_k1], h2_asymp_se[last_k2]), 3),
  ESS_pred        = round(c(h1_ess_pred[last_k1], h2_ess_pred[last_k2]), 1),
  Pred_Eff_Ratio  = round(c(h1_ess_pred[last_k1], h2_ess_pred[last_k2]) / S_last, 4),
  min_ESS_param   = round(c(h1_min_ess[last_k1], h2_min_ess[last_k2]), 2),
  Param_Eff_Ratio = format(c(h1_min_ess[last_k1], h2_min_ess[last_k2]) / S_last,
                           scientific = TRUE, digits = 3)
)

cat("--- Model Convergence Summary at S_max ---\n")
print(summary_df, row.names = FALSE)
cat("\n")

delta_elpd <- h1_elpd[last_k1] - h2_elpd[last_k2]
mcse_delta <- sqrt(h1_asymp_se[last_k1]^2 + h2_asymp_se[last_k2]^2)
z_delta    <- abs(delta_elpd) / mcse_delta
p_delta    <- 2 * (1 - pnorm(z_delta))
ci_lo      <- delta_elpd - 1.96 * mcse_delta
ci_hi      <- delta_elpd + 1.96 * mcse_delta

cat("--- Out-of-Sample Predictive Comparison ---\n")
cat(sprintf("  Delta ELPD (M1 - M2) : %+10.3f\n", delta_elpd))
cat(sprintf("  Monte Carlo SE       : %10.3f\n", mcse_delta))
cat(sprintf("  95%% MC Confidence Int: [%.3f, %.3f]\n", ci_lo, ci_hi))
cat(sprintf("  Z-score              : %10.2f (p = %s)\n", z_delta, format.pval(p_delta, digits = 3)))
if (delta_elpd > 0 && ci_lo > 0) {
  cat("  Conclusion           : Model M1 (Joint) significantly outperforms Model M2 (Separate)!\n\n")
} else if (delta_elpd < 0 && ci_hi < 0) {
  cat("  Conclusion           : Model M2 (Separate) significantly outperforms Model M1 (Joint)!\n\n")
} else {
  cat("  Conclusion           : No statistically significant difference detected.\n\n")
}

## 95% Monte Carlo intervals
z_crit <- qnorm(0.975)

m1_lo <- h1_elpd - z_crit * h1_asymp_se
m1_hi <- h1_elpd + z_crit * h1_asymp_se
m2_lo <- h2_elpd - z_crit * h2_asymp_se
m2_hi <- h2_elpd + z_crit * h2_asymp_se

## Plot styling and axes
COL_M1     <- "#08519C"
COL_M2     <- "#D73027"
COL_GRID   <- "grey92"
COL_BOX    <- "grey75"

x_lo        <- if (from_zero) 0 else S_star
xlim_common <- c(x_lo, S_max)
s_ticks     <- pretty(xlim_common, n = 6)
s_ticks     <- s_ticks[s_ticks >= x_lo & s_ticks <= S_max]
if (length(s_ticks) < 2) s_ticks <- S_vals
s_labels    <- if (S_max >= 1000) ifelse(s_ticks == 0, "0", paste0(s_ticks / 1000, "K")) else as.character(s_ticks)

elpd_all_lo  <- min(c(m1_lo, m2_lo), na.rm = TRUE)
elpd_all_hi  <- max(c(m1_hi, m2_hi), na.rm = TRUE)
y_range_elpd <- elpd_all_hi - elpd_all_lo
if (y_range_elpd < 1e-4) y_range_elpd <- 10.0

## y-limits hug the min / max of the plotted intervals (extended only if the legend needs it)
ylim_tight <- c(elpd_all_lo, elpd_all_hi) + c(-1, 1) * 0.04 * y_range_elpd

## The legend goes in the top-right or bottom-right corner of the panel. Starting from
## y-limits that hug the data, the corner already clear of the curves is used; if the
## curves reach into both, the y-axis is extended on one side by just enough for the
## legend to sit clear of them. Returns the corner and the final y-limits.
LEG_INSET <- c(0.015, 0.025)
fit_legend_corner <- function(xlim, ylim, x_list, lo_list, hi_list, leg_args, gap = 0.03) {
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, yaxs = "i")
  usr <- par("usr")
  r   <- do.call(legend, c(list("topright", inset = LEG_INSET, plot = FALSE), leg_args))$rect
  fh  <- (usr[4] - (r$top - r$h)) / (usr[4] - usr[3]) + gap   # share of panel height needed
  under <- unlist(x_list) >= r$left - 0.02 * (usr[2] - usr[1]) # checkpoints below / above the box
  d_hi  <- max(unlist(hi_list)[under], na.rm = TRUE)
  d_lo  <- min(unlist(lo_list)[under], na.rm = TRUE)
  top_ext <- max(0, (d_hi - fh * ylim[1]) / (1 - fh) - ylim[2])
  bot_ext <- max(0, ylim[1] - (d_lo - fh * ylim[2]) / (1 - fh))
  par(new = TRUE)                                              # the real plot reuses this page
  if (top_ext <= bot_ext) list(pos = "topright",    ylim = c(ylim[1], ylim[2] + top_ext))
  else                    list(pos = "bottomright", ylim = c(ylim[1] - bot_ext, ylim[2]))
}

## ELJPD tracking with 95% intervals
render_elpd_panel <- function(legend_cex = 1.40) {
  leg_args <- list(
         legend = c(
           expression(italic(M)[1] ~ "(Joint Covariance GP)"),
           expression(italic(M)[2] ~ "(Separate Covariance GP)")
         ),
         col    = c(COL_M1, COL_M2),
         pch    = c(19, 1),
         lty    = c(1, 1),
         lwd    = c(3.2, 3.2),
         pt.cex = c(1.7, 1.7),
         pt.lwd = c(NA, 2.2),
         bty    = "o", bg = adjustcolor("white", 0.94), box.col = COL_BOX,
         cex    = legend_cex, y.intersp = 1.15)
  fit <- fit_legend_corner(xlim_common, ylim_tight, list(S1_vals, S2_vals),
                           list(m1_lo, m2_lo), list(m1_hi, m2_hi), leg_args)
  ylim_elpd     <- fit$ylim
  y_ticks_elpd  <- pretty(ylim_elpd, n = 6)
  y_ticks_elpd  <- y_ticks_elpd[y_ticks_elpd >= ylim_elpd[1] & y_ticks_elpd <= ylim_elpd[2]]
  y_labels_elpd <- format(y_ticks_elpd, big.mark = ",")

  plot(NA, xlim = xlim_common, ylim = ylim_elpd, yaxs = "i",
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = "", main = "", xaxt = "n", yaxt = "n")
  title(ylab = expression(widehat(ELJPD)(italic(S))), line = 8)

  axis(1, at = s_ticks, labels = s_labels, gap.axis = 0)
  axis(2, at = y_ticks_elpd, labels = y_labels_elpd)
  abline(h = y_ticks_elpd, col = COL_GRID, lty = 1, lwd = 0.8)
  abline(v = s_ticks,      col = COL_GRID, lty = 1, lwd = 0.8)

  arrows(S1_vals, m1_lo, S1_vals, m1_hi,
         code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M1, 0.45), lwd = 2.4)
  arrows(S2_vals, m2_lo, S2_vals, m2_hi,
         code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M2, 0.50), lwd = 2.4)

  lines(S1_vals,  h1_elpd, col = COL_M1, lwd = 3.5, lty = 1)
  points(S1_vals, h1_elpd, col = COL_M1, pch = 19, cex = 1.7)
  lines(S2_vals,  h2_elpd, col = COL_M2, lwd = 3.5, lty = 1)
  points(S2_vals, h2_elpd, col = COL_M2, pch = 1, cex = 1.7, lwd = 2.2)

  do.call(legend, c(list(fit$pos, inset = LEG_INSET), leg_args))
  box(lwd = 1.8)
}

## Export (same size as the arrhythmia figure, output_arrhythmia.R)
FIG_WIDTH  <- 10
FIG_HEIGHT <- 7.5
fig_file <- paste0("elpd_mcse_spatial", run_tag, ".pdf")
pdf(fig_file, width = FIG_WIDTH, height = FIG_HEIGHT, useDingbats = FALSE)
par(
  mar      = c(6.5, 11, 1.2, 1.5),
  oma      = c(0.5, 0.5, 0.5, 0.5),
  mgp      = c(5.2, 1.5, 0),
  cex.lab  = 2.6,
  cex.axis = 2.1,
  las      = 1
)
render_elpd_panel(legend_cex = 1.45)
invisible(dev.off())
cat(sprintf("Saved: %s\n", fig_file))
