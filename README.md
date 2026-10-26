# Monte Carlo uncertainty quantification for Bayesian cross-validation

Code and saved outputs for the simulation studies and real-data analyses in the paper
*"Monte Carlo uncertainty quantification for Bayesian cross-validation"*.

In every study, two competing models (M1 and M2) are compared by cross-validated predictive
accuracy. Each model is fitted with a single MCMC chain that targets the mixture proposal
`q_mix`, and the CV estimate (ELPPD for LOO-CV, ELJPD for K-fold CV) is computed by mixture
importance sampling. The asymptotic Monte Carlo standard error (MCSE) and the effective
sample size of the estimator are computed with batch means (`mcmcse`). They are used to
(i) track the estimate and its 95% Monte Carlo interval as the chain length `S` grows, together
with the separation point `S_sep` after which the intervals of M1 and M2 no longer overlap,
and (ii) stop the chain once a target precision `eps` is reached.

## Repository structure

```
Simulation_studies/                        Section 6
├── Bayesian_linear_regression/            Section 6.1  Subset selection in Bayesian linear regression
│   ├── non_replicated_study/              single chain per model
│   ├── replicated_study/                  500 independent chains per model
│   └── output_analysis/                   figure: elpd_mcse_linreg.pdf
└── Smoothness_misspec_GP/                 Section 6.2  Smoothness misspecification in GP regression
    ├── non_replicated_study/              single chain per model, exact and brute-force references
    ├── replicated_study/                  500 independent chains per model
    └── output_analysis/                   figure: elpd_mcse_smoothness.pdf

Real_data_analysis/                        Section 7
├── Arrhythmia_LOOCV/                      Section 7.1  Link misspecification in binary regression
│                                          figure: elpd_mcse_arrhythmia.pdf
└── Mixed_type_GP/                         Section 7.2  Latent covariance misspecification in mixed-type spatial data
    ├── M1_estimation/                     joint (cross-covariance) GP model
    ├── M2_estimation/                     independent GP model
                                           figure: elpd_mcse_spatial.pdf
```

Each study folder contains the scripts that run the samplers (`run_*.R`, with C++ samplers in
`*.cpp`), the saved sampler outputs (`*_outputs.RData`, `*_outputs_reps.RData`), and a
plotting script that reads those outputs and writes the figures. The `.Rout` files are the R
console logs of the runs that produced the saved real-data and replicated GP outputs.

## Studies

### Section 6.1: Subset selection in Bayesian linear regression

`n = 500`, `d = 250`, Gaussian prior on the coefficients. M1 is the true model and M2 drops one
active predictor; the coefficient of the dropped predictor is calibrated so that the exact
LOO-ELPD gap equals 4 (`sim_data.R`). The exact LOO-ELPD of both models is available in closed
form and is used as a reference.

| File | Purpose |
| --- | --- |
| `*/sim_data.R` | Simulates the data and computes the exact LOO-ELPD |
| `*/run_true_d250_n500.R`, `*/run_misspec_d250_n500.R` | Sampler runs for M1 and M2 |
| `*/elppd_utils.R`, `proposals.R`, `sampler.R` | Mixture IS estimator, mixture target, sampler |
| `output_analysis/output_analysis.R` | Produces `elpd_mcse_linreg.pdf` |

### Section 6.2: Smoothness misspecification in GP regression

Vecchia-approximated GP regression with a nugget on a 20 x 20 grid, 10-fold CV. M1 uses the true
Matérn smoothness `nu = 0.5`; M2 uses `nu = 1`.

| File | Purpose |
| --- | --- |
| `*/sim_data.R` | Simulates the data (`sim_data.RData`) |
| `*/gibbs_kfold.cpp` | Gibbs sampler on `q_mix` |
| `*/model_utils.R`, `vecchia.R`, `eljpd.utils.R` | Set-up, Vecchia neighbour sets, ELJPD and MCSE |
| `non_replicated_study/run_M1.R`, `run_M2.R` | Single-chain runs (`M1_outputs.RData`, `M2_outputs.RData`) |
| `non_replicated_study/exact_reference.R` | Exact 10-fold ELJPD by quadrature (`exact_reference.RData`) |
| `non_replicated_study/run_bruteforce.R` | Optional brute-force reference (one chain per fold) |
| `replicated_study/run_M1.R`, `run_M2.R` | 500 replicated chains (`M1_outputs_reps.RData`, `M2_outputs_reps.RData`) |
| `output_analysis/output_analysis.R` | Produces `elpd_mcse_smoothness.pdf` |

### Section 7.1: Link misspecification in binary regression

Arrhythmia data (`n = 452`, `d = 258` including the intercept) with a regularised horseshoe prior.
M1 uses the logistic link and M2 the complementary log-log link. Model comparison is by LOO-CV.

| File | Purpose |
| --- | --- |
| `data.RData`, `model_specs.RData` | Design matrix and binary response |
| `gibbs_loocv_ess.cpp` | Elliptical-slice-within-Gibbs sampler on `q_mix` |
| `elppd_utils.R` | Mixture IS estimator of the ELPPD and its MCSE |
| `run_M1.R`, `run_M2.R` | Sampler runs (`M1_outputs.RData`, `M2_outputs.RData`) |
| `output_arrhythmia.R` | Produces `elpd_mcse_arrhythmia.pdf` |

### Section 7.2: Latent covariance misspecification in mixed-type spatial data

US wildfire data on `N = 3,503` grid cells with a Gaussian response (log burned area) and a Poisson
response (fire counts). M1 links the two latent processes with a separable cross-covariance GP;
M2 uses independent GPs. Both use a Vecchia approximation, and comparison is by 10-fold CV.

| File | Purpose |
| --- | --- |
| `data_full.RData` | Wildfire data |
| `M1_estimation/gibbs_kfold_joint.cpp`, `M2_estimation/gibbs_kfold_sep.cpp` | Gibbs samplers on `q_mix` |
| `M*_estimation/vecchia.R`, `eljpd.utils.R` | Max-min ordering, neighbour sets, ELJPD and MCSE |
| `M1_estimation/run_M1.R`, `M2_estimation/run_M2.R` | Sampler runs (`M1_outputs.RData`, `M2_outputs.RData`) |
| `output_mixed_type_spatial.R` | Produces `elpd_mcse_spatial.pdf` |

## Reproducing the figures

The saved outputs are included, so each figure can be regenerated in a few seconds from the
plotting script in its folder:

```r
source("Simulation_studies/Bayesian_linear_regression/output_analysis/output_analysis.R")
source("Simulation_studies/Smoothness_misspec_GP/output_analysis/output_analysis.R")
source("Real_data_analysis/Arrhythmia_LOOCV/output_arrhythmia.R")
source("Real_data_analysis/Mixed_type_GP/output_mixed_type_spatial.R")
```

Each script sets the working directory to its own folder (via `this.path`) and writes its
figures there. All figures are exported as PDF only.

To regenerate the sampler outputs, run the `run_*.R` scripts in each folder before the
plotting script. The replicated studies run chains in parallel (`foreach`/`doParallel`; the
GP scripts use a FORK cluster) and are computationally intensive. Set `n_cores` to suit the
machine.

## Requirements

R (≥ 4.5) with a C++ toolchain for `Rcpp`, and the packages

```r
install.packages(c("Rcpp", "RcppArmadillo", "mcmcse", "coda", "fields", "FNN", "sp", "gstat",
                   "glmnet", "MASS", "foreach", "doParallel", "this.path"))
```
