# ==============================================================================
# Reproducible R Script to Generate Publication-Ready LaTeX Table for ELPPD Variance
# ==============================================================================

# Set working directory to the script's directory (if run from a file/RStudio)
if (requireNamespace("this.path", quietly = TRUE)) {
  setwd(this.path::here())
}

# Configuration
K.vec <- c(5, 10, 50)
S_mcmc <- 5e5 # MCMC Sample size (iterations)

# Methods mapping to their publication labels
methods_labels <- c(
  g_is   = "Gelfand IS",
  g_tis  = "Gelfand TIS",
  g_psis = "Gelfand PSIS",
  b_is   = "Bronze IS",
  b_tis  = "Bronze TIS",
  b_psis = "Bronze PSIS",
  z_is   = "Zanella IS",
  z_tis  = "Zanella TIS",
  z_psis = "Zanella PSIS"
)

# Start building the LaTeX tabular structure
latex_code <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Comparison of Sample Variance and Estimated Variance (on the Asymptotic Scale $S \\times \\text{Variance}$) of the $K$-fold ELPPD Estimator. All variance values are scaled by the MCMC sample size $S = 500,000$.}",
  "\\label{tab:elppd_variance_comparison}",
  "\\begin{tabular}{clcccc}",
  "\\toprule",
  "\\multirow{2.5}{*}{$K$} & \\multirow{2.5}{*}{Method} & Sample Variance & \\multicolumn{2}{c}{Estimated Variance ($S \\times \\widehat{\\text{Var}}_{\\text{est}}$)} & Ratio \\\\",
  "\\cmidrule(lr){4-5}",
  " & & ($S \\times \\text{Var}(\\widehat{\\text{ELPPD}})$) & Median & 95\\% CI & (Sample / Est) \\\\",
  "\\midrule"
)

# Loop through each fold setting and method
for (K in K.vec) {
  file_name <- paste0(K, "_fold_elppd_Gaussian.RData")
  
  if (!file.exists(file_name)) {
    stop(paste("File", file_name, "not found in current directory:", getwd()))
  }
  
  load(file_name) # Loads 'results' list
  
  first_in_block <- TRUE
  
  for (m in names(methods_labels)) {
    var_col <- paste0("var_", m)
    elppd_col <- paste0("elppd_", m)
    
    # Scale both estimated and sample variances to the asymptotic scale (multiplied by S)
    est_vars <- results[[var_col]] * S_mcmc
    samp_var <- var(results[[elppd_col]]) * S_mcmc
    
    median_est <- median(est_vars)
    ci_est <- quantile(est_vars, probs = c(0.025, 0.975))
    ratio <- samp_var / median_est
    
    # Format the multirow label for K on the first row of each block
    k_label <- if (first_in_block) paste0("\\multirow{9.5}{*}{", K, "}") else ""
    first_in_block <- FALSE
    
    # Clean numerical formatting for journal output and wrap in LaTeX math mode ($)
    fmt <- function(x) {
      if (is.na(x)) return("NA")
      if (x > 1e4 || x < 1e-2) {
        # Format using LaTeX scientific notation syntax: $a \times 10^{b}$
        s <- formatC(x, format = "e", digits = 2)
        parts <- strsplit(s, "e")[[1]]
        exp_val <- as.numeric(parts[2])
        paste0("$", parts[1], " \\times 10^{", exp_val, "}$")
      } else {
        paste0("$", formatC(x, format = "f", digits = 2), "$")
      }
    }
    
    ci_str <- paste0("[", fmt(ci_est[1]), ", ", fmt(ci_est[2]), "]")
    
    # Generate the LaTeX row
    row_str <- paste(
      k_label, "&",
      methods_labels[m], "&",
      fmt(samp_var), "&",
      fmt(median_est), "&",
      ci_str, "&",
      formatC(ratio, format = "f", digits = 2),
      "\\\\"
    )
    
    latex_code <- c(latex_code, row_str)
  }
  
  # Add midrule separator between groups (except the last group)
  if (K != tail(K.vec, 1)) {
    latex_code <- c(latex_code, "\\midrule")
  }
}

# Close the LaTeX code
latex_code <- c(
  latex_code,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 1. Save table directly to file
writeLines(latex_code, "variance_comparison_table.tex")
message("Saved table to 'variance_comparison_table.tex'")

# 2. Output code directly to console for copy-pasting
cat("\n--- Generated LaTeX Table Code ---\n\n")
cat(paste(latex_code, collapse = "\n"), "\n")