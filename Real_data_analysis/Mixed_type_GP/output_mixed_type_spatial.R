## Latent covariance misspecification in mixed-type spatial data (Section 7.2): US wildfires, N = 3,503 grid cells
## M1: joint (separable cross-covariance) GP, M2: independent GPs; output elpd_mcse_spatial.pdf

rm(list = ls())
invisible(gc())

if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

## Load model outputs
cat("==============================================================================\n")
cat("Real Data Analysis: Model Comparison & Convergence Diagnostics\n")
cat("==============================================================================\n\n")

file_m1_candidates <- c("M1_outputs.RData", file.path("M1_estimation", "M1_outputs.RData"))
file_m2_candidates <- c("M2_outputs.RData", file.path("M2_estimation", "M2_outputs.RData"))

file_m1 <- file_m1_candidates[file.exists(file_m1_candidates)][1]
file_m2 <- file_m2_candidates[file.exists(file_m2_candidates)][1]

if (is.na(file_m1) || is.na(file_m2)) {
  stop("Output files missing! Please ensure M1_outputs.RData and M2_outputs.RData are present.")
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

h1_elpd     <- extract_field(h1, c("eljpd_hat", "elpd_hat"))
h1_asymp_se <- extract_field(h1, c("asymp_se", "se"))
h1_ess_pred <- extract_field(h1, c("ess_elpd", "ess"))
h1_min_ess  <- extract_field(h1, "minESS_xi", default = rep(NA_real_, nrow(h1)))

h2_elpd     <- extract_field(h2, c("eljpd_hat", "elpd_hat"))
h2_asymp_se <- extract_field(h2, c("asymp_se", "se"))
h2_ess_pred <- extract_field(h2, c("ess_elpd", "ess"))
h2_min_ess  <- extract_field(h2, "minESS_xi", default = rep(NA_real_, nrow(h2)))

S_vals        <- h1$S
S_star        <- get_var(env_m1, "S_star", min(S_vals))
S_max         <- get_var(env_m1, "S_max", max(S_vals))
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
cat(sprintf("Metropolis acceptance rate (phi): M1 = %.4f, M2 = %.4f\n\n", acc_m1, acc_m2))

## Summary at S_max
last_k <- length(S_vals)

summary_df <- data.frame(
  Model           = c("M1 (Joint)", "M2 (Separate)"),
  Specification   = c("Separable Cross-Covariance GP", "Independent Univariate GPs"),
  S_max           = rep(format(as.integer(S_vals[last_k]), big.mark = ","), 2),
  ELPD_hat        = round(c(h1_elpd[last_k], h2_elpd[last_k]), 2),
  MCSE            = round(c(h1_asymp_se[last_k], h2_asymp_se[last_k]), 3),
  ESS_pred        = round(c(h1_ess_pred[last_k], h2_ess_pred[last_k]), 1),
  Pred_Eff_Ratio  = round(c(h1_ess_pred[last_k], h2_ess_pred[last_k]) / S_vals[last_k], 4),
  min_ESS_param   = round(c(h1_min_ess[last_k], h2_min_ess[last_k]), 2),
  Param_Eff_Ratio = format(c(h1_min_ess[last_k], h2_min_ess[last_k]) / S_vals[last_k],
                           scientific = TRUE, digits = 3)
)

cat("--- Model Convergence Summary at S_max ---\n")
print(summary_df, row.names = FALSE)
cat("\n")

delta_elpd <- h1_elpd[last_k] - h2_elpd[last_k]
mcse_delta <- sqrt(h1_asymp_se[last_k]^2 + h2_asymp_se[last_k]^2)
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

## S_sep: earliest S after which the model leading at S_max keeps a disjoint 95% interval
z_crit <- qnorm(0.975)

m1_lo <- h1_elpd - z_crit * h1_asymp_se
m1_hi <- h1_elpd + z_crit * h1_asymp_se
m2_lo <- h2_elpd - z_crit * h2_asymp_se
m2_hi <- h2_elpd + z_crit * h2_asymp_se

m1_wins <- h1_elpd[last_k] >= h2_elpd[last_k]

compute_S_sep <- function(S, lo1, hi1, lo2, hi2, m1_wins) {
  is_sep <- if (m1_wins) lo1 > hi2 else lo2 > hi1
  if (length(is_sep) > 0 && is_sep[length(is_sep)]) {
    fi <- which(!is_sep)
    k_sep <- if (length(fi) == 0) 1L else max(fi) + 1L
    return(S[k_sep])
  }
  NA_real_
}

sep_S <- compute_S_sep(S_vals, m1_lo, m1_hi, m2_lo, m2_hi, m1_wins)

winner_name <- if (m1_wins) "M1 (Joint)" else "M2 (Separate)"
loser_name  <- if (m1_wins) "M2 (Separate)" else "M1 (Joint)"
winner_lo   <- if (m1_wins) m1_lo else m2_lo
loser_hi    <- if (m1_wins) m2_hi else m1_hi

cat("--- Empirical Separation Point Analysis ---\n")
if (!is.na(sep_S)) {
  k_sep     <- which(S_vals == sep_S)
  gap_sep   <- winner_lo[k_sep] - loser_hi[k_sep]
  sep_k_str <- if (sep_S >= 1000) sprintf(" (%sK)", format(as.integer(sep_S / 1000))) else ""
  cat(sprintf("  Superior model       : %s\n", winner_name))
  cat(sprintf("  Empirical separation attained at S_sep = %s iterations%s!\n",
              format(as.integer(sep_S), big.mark = ","), sep_k_str))
  cat(sprintf("  At S = %s:\n", format(as.integer(sep_S), big.mark = ",")))
  cat(sprintf("    %s lower 95%% bound : %.2f\n", winner_name, winner_lo[k_sep]))
  cat(sprintf("    %s upper 95%% bound : %.2f\n", loser_name, loser_hi[k_sep]))
  cat(sprintf("    Strict separation gap : +%.2f (ELPD units)\n", gap_sep))
  cat(sprintf("  Separation holds permanently from S = %s through S_max = %s.\n",
              format(as.integer(sep_S), big.mark = ","),
              format(as.integer(S_max), big.mark = ",")))
  cat("\n")
} else {
  cat("  Confidence intervals continue to overlap up to S_max; no complete separation attained.\n\n")
}

## Plot styling and axes
COL_M1     <- "#08519C"
COL_M2     <- "#D73027"
COL_SEP    <- "#006837"
COL_GRID   <- "grey92"
COL_BOX    <- "grey75"

xlim_common <- range(S_vals)
s_ticks     <- pretty(S_vals, n = 6)
s_ticks     <- s_ticks[s_ticks >= S_star & s_ticks <= S_max]
if (length(s_ticks) < 2) s_ticks <- S_vals
s_labels    <- if (S_max >= 1000) paste0(s_ticks / 1000, "K") else as.character(s_ticks)

elpd_all_lo  <- min(c(m1_lo, m2_lo), na.rm = TRUE)
elpd_all_hi  <- max(c(m1_hi, m2_hi), na.rm = TRUE)
y_range_elpd <- elpd_all_hi - elpd_all_lo
if (y_range_elpd < 1e-4) y_range_elpd <- 10.0

y_ticks_elpd  <- pretty(c(elpd_all_lo - 0.05 * y_range_elpd, elpd_all_hi + 0.32 * y_range_elpd), n = 7)
ylim_elpd     <- range(y_ticks_elpd)
y_labels_elpd <- format(y_ticks_elpd, big.mark = ",")

## ELJPD tracking with 95% intervals and separation point
render_elpd_panel <- function(legend_cex = 1.40) {
  plot(NA, xlim = xlim_common, ylim = ylim_elpd,
       xlab = expression("Iterations (" * italic(S) * ")"),
       ylab = "", main = "", xaxt = "n", yaxt = "n")
  title(ylab = expression(widehat(ELJPD)(italic(S))), line = 8)

  axis(1, at = s_ticks, labels = s_labels, gap.axis = 0)
  axis(2, at = y_ticks_elpd, labels = y_labels_elpd)
  abline(h = y_ticks_elpd, col = COL_GRID, lty = 1, lwd = 0.8)
  abline(v = s_ticks,      col = COL_GRID, lty = 1, lwd = 0.8)

  if (!is.na(sep_S)) abline(v = sep_S, col = COL_SEP, lty = 4, lwd = 3.2)

  arrows(S_vals, m1_lo, S_vals, pmin(m1_hi, ylim_elpd[2]),
         code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M1, 0.45), lwd = 2.4)
  arrows(S_vals, m2_lo, S_vals, pmin(m2_hi, ylim_elpd[2]),
         code = 3, angle = 90, length = 0.06,
         col = adjustcolor(COL_M2, 0.50), lwd = 2.4)

  lines(S_vals,  h1_elpd, col = COL_M1, lwd = 3.5, lty = 1)
  points(S_vals, h1_elpd, col = COL_M1, pch = 19, cex = 1.7)
  lines(S_vals,  h2_elpd, col = COL_M2, lwd = 3.5, lty = 1)
  points(S_vals, h2_elpd, col = COL_M2, pch = 1, cex = 1.7, lwd = 2.2)

  sep_txt <- if (!is.na(sep_S)) {
    bquote("Separation (" * italic(S)[sep] == .(format(as.integer(sep_S), big.mark = ",")) * ")")
  } else {
    expression("No separation")
  }

  legend("topright", inset = c(0.025, 0.025),
         legend = c(
           expression(italic(M)[1] ~ "(Joint Covariance GP)"),
           expression(italic(M)[2] ~ "(Separate Covariance GP)"),
           sep_txt
         ),
         col    = c(COL_M1, COL_M2, COL_SEP),
         pch    = c(19, 1, NA),
         lty    = c(1, 1, 4),
         lwd    = c(3.2, 3.2, 3.2),
         pt.cex = c(1.7, 1.7, NA),
         pt.lwd = c(NA, 2.2, NA),
         bty    = "o", bg = adjustcolor("white", 0.94), box.col = COL_BOX,
         cex    = legend_cex, y.intersp = 1.15)
  box(lwd = 1.8)
}

## Export
pdf("elpd_mcse_spatial.pdf", width = 10, height = 9.5, useDingbats = FALSE)
par(
  mar      = c(6.5, 11, 2.5, 2.0),
  oma      = c(0.5, 0.5, 0.5, 0.5),
  mgp      = c(5.2, 1.5, 0),
  cex.lab  = 2.6,
  cex.axis = 2.1,
  las      = 1
)
render_elpd_panel(legend_cex = 1.45)
invisible(dev.off())
cat("Saved: elpd_mcse_spatial.pdf\n")
