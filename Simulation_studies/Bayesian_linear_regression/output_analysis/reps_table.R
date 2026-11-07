## Subset selection in Bayesian linear regression (Section 6.1): replicated study summary table
## For each replication r and precision eps, the chain is terminated at the first checkpoint
## S_term(r) with ESS_ELPPD >= 4 z^2 / eps^2 (S_max if never reached). All metrics, including
## coverage of the nominal 95% asymptotic interval, are evaluated at that replication's S_term.
## Output: reps_summary_linreg.tex (LaTeX table via xtable)

rm(list = ls())

library(xtable)

if (requireNamespace("this.path", quietly = TRUE)) setwd(this.path::here())

rep_dir  <- "../replicated_study"
tex_file <- "reps_summary_linreg.tex"

alpha    <- 0.05
z_crit   <- qnorm(1 - alpha / 2)
eps_vals <- c(0.05, 0.025)
ess_tgt  <- 4 * z_crit^2 / eps_vals^2

load_env <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e
}

e1 <- load_env(file.path(rep_dir, "M1_outputs_reps.RData"))
e2 <- load_env(file.path(rep_dir, "M2_outputs_reps.RData"))

reps_list <- list(M1 = e1$results_reps, M2 = e2$results_reps)
ref       <- c(M1 = e1$elpd_M1, M2 = e1$elpd_M2)   # exact analytical LOO-ELPD
R         <- length(reps_list$M1)

# Termination point of one replication; asymp_se is stored as the MC variance
get_term <- function(h, target) {
  ess <- if (!is.null(h$ess_elppd)) h$ess_elppd else h$ess
  k   <- which(ess >= target)[1]
  reached <- !is.na(k)
  if (!reached) k <- length(h$S)
  data.frame(S_term   = h$S[k],
             ELPPD    = h$elppd_hat[k],
             half_len = z_crit * sqrt(max(h$asymp_se[k], 0)),
             reached  = reached)
}

summ <- do.call(rbind, lapply(names(reps_list), function(mod) {
  do.call(rbind, lapply(seq_along(eps_vals), function(j) {
    d  <- do.call(rbind, lapply(reps_list[[mod]], function(x) get_term(x$history, ess_tgt[j])))
    cov <- abs(d$ELPPD - ref[mod]) <= d$half_len
    q   <- quantile(d$S_term, c(0.25, 0.5, 0.75), names = FALSE)
    data.frame(model = mod, eps = eps_vals[j],
               S_med = q[2], S_q1 = q[1], S_q3 = q[3],
               half_len = mean(d$half_len),
               coverage = 100 * mean(cov),
               not_reached = sum(!d$reached))
  }))
}))

cat(sprintf("Exact LOO-ELPD  M1: %.2f  M2: %.2f  difference: %.2f  (R = %d)\n",
            ref["M1"], ref["M2"], ref["M1"] - ref["M2"], R))
print(summ, digits = 4, row.names = FALSE)
if (any(summ$not_reached > 0))
  warning("Some replications never reached the ESS target; S_max was used for those.")

# LaTeX formatting helpers
fmt_int <- function(x) formatC(x, format = "d", big.mark = "{,}")
fmt_K   <- function(x) paste0(sprintf("%g", x / 1000), "\\text{K}")

tab <- data.frame(
  Model = ifelse(duplicated(summ$model), "",
                 sprintf("\\multirow{2}{*}{$\\mc{M}_%s$}", sub("M", "", summ$model))),
  Eps   = sprintf("$%.3f$", summ$eps),
  Sterm = sprintf("$%s$ [$%s, %s$]", fmt_int(round(summ$S_med)),
                  fmt_K(summ$S_q1), fmt_K(summ$S_q3)),
  Half  = sprintf("$%.3f$", summ$half_len),
  Cov   = sprintf("$%.1f\\%%$", summ$coverage),
  stringsAsFactors = FALSE
)
colnames(tab) <- c("Model", "Precision $\\varepsilon$", "Median $\\Sterm$ [IQR]",
                   "$\\text{half-len}(\\bar{C_\\alpha})$", "Coverage")

caption <- sprintf(paste0(
  "Operational termination metrics and empirical coverage of nominal %d\\%% asymptotic ",
  "confidence intervals across $R = %d$ independent MCMC replications. The analytical ",
  "reference values are $\\mathrm{ELPD}_{\\mc{D}}(\\mc{M}_1) = %.2f$ and ",
  "$\\mathrm{ELPD}_{\\mc{D}}(\\mc{M}_2) = %.2f$, giving $\\Delta\\mathrm{ELPD}_{\\mc{D}} = %s$."),
  round(100 * (1 - alpha)), R, ref["M1"], ref["M2"],
  sprintf("%g", round(ref["M1"] - ref["M2"], 2)))

xt <- xtable(tab, caption = caption, label = "tab:reps_summary", align = "llcccc")

# \midrule between models (after the M1 rows)
n_M1 <- sum(summ$model == "M1")
print(xt,
      file = tex_file,
      table.placement = "H",
      caption.placement = "top",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      size = "\\small\\setlength{\\tabcolsep}{5.5pt}",
      add.to.row = list(pos = list(n_M1), command = "\\midrule\n"),
      comment = FALSE)

cat(readLines(tex_file), sep = "\n")
cat(sprintf("\nSaved %s\n", tex_file))
