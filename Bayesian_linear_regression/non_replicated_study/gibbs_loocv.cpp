#include <Rcpp.h>
using namespace Rcpp;

// =============================================================================
// HMC-within-Gibbs sampler for the LOOCV mixture proposal
//
//     q_mix(theta, k)  ∝  p(theta) p(y_{-k} | theta),      k = 1, ..., n
//
// One Gibbs iteration:
//   1. k | theta  ~  Categorical( p(y_k | theta)^{-1} / sum_j p(y_j | theta)^{-1} )
//        drawn on the LOG scale with the Gumbel-max trick:
//            k = argmax_i { -log p(y_i | theta) + G_i },   G_i ~ Gumbel(0, 1)
//        (no exponentials, no normalising constant, no under/overflow).
//   2. theta | k  ~  p(theta | y_{-k})  : n_hmc HMC transitions that leave
//        p(theta | y_{-k}) invariant (Metropolis-within-Gibbs), started at the
//        current theta.
//
// Discarding k, the theta draws are draws from the marginal q_mix(theta).
// Same arguments / return names as hmc_loocv_cpp() (plus n_hmc and the fold labels).
// =============================================================================

// Numerically stable log-sum-exp
// [[Rcpp::export]]
double log_sum_exp(const NumericVector& x) {
  int n = x.size();
  if (n == 0) return R_NegInf;
  double m = x[0];
  for (int i = 1; i < n; ++i) {
    if (x[i] > m) m = x[i];
  }
  if (R_IsNA(m) || R_IsNaN(m)) return NA_REAL;
  if (m == R_PosInf) return R_PosInf;
  if (m == R_NegInf) return R_NegInf;
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += std::exp(x[i] - m);
  }
  return m + std::log(s);
}


// resid = y - X theta      (X is column-major: contiguous column access)
static void residuals(const double* X, const double* y, const double* theta,
                      int n, int d, double* resid) {
  for (int i = 0; i < n; ++i) resid[i] = y[i];
  for (int j = 0; j < d; ++j) {
    const double* col = X + (std::size_t)j * n;
    const double t = theta[j];
    for (int i = 0; i < n; ++i) resid[i] -= col[i] * t;
  }
}

// Potential of the fold-k training posterior (up to an additive constant)
//     U_k(theta) = -log p(theta) - sum_{i != k} log p(y_i | theta)
// and its gradient. `resid` is a scratch vector of length n.
static double potential_fold(const double* X, const double* y, const double* theta,
                             int n, int d, double sigma2, double tau2, int k,
                             double* resid, double* grad) {
  residuals(X, y, theta, n, d, resid);
  resid[k] = 0.0;                                  // leave observation k out

  double ssq = 0.0;
  for (int i = 0; i < n; ++i) ssq += resid[i] * resid[i];

  double th2 = 0.0;
  const double inv_s = 1.0 / sigma2;
  for (int j = 0; j < d; ++j) {
    const double* col = X + (std::size_t)j * n;
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += col[i] * resid[i];
    grad[j] = theta[j] / tau2 - s * inv_s;         // grad of U_k
    th2 += theta[j] * theta[j];
  }
  return 0.5 * (th2 / tau2 + ssq * inv_s);
}

// [[Rcpp::export]]
List gibbs_loocv_cpp(
    const NumericMatrix& X,
    const NumericVector& y,
    const NumericVector& theta_init,
    double sigma2,
    double tau2,
    double epsilon,
    int L,
    int niters,
    int n_hmc = 1            // HMC transitions per Gibbs iteration
) {
  const int n = X.nrow();
  const int d = X.ncol();
  const double* Xp = &X(0, 0);
  const double* yp = &y[0];
  const double inv_2s = 0.5 / sigma2;

  NumericMatrix samples(niters, d);
  IntegerVector folds(niters);                     // sampled fold labels (1-based)
  std::vector<double> theta(d), theta_p(d), p(d), p_p(d), g_cur(d), g(d);
  std::vector<double> resid(n), score(n);

  for (int j = 0; j < d; ++j) theta[j] = theta_init[j];
  double acc = 0.0;

  for (int iter = 0; iter < niters; ++iter) {

    // ---- Step 1: k | theta  (Gumbel-max on the log scale) -------------------
    // -log p(y_i | theta) = const + r_i^2 / (2 sigma2); the constant cancels.
    residuals(Xp, yp, theta.data(), n, d, resid.data());
    int k = 0;
    double best = R_NegInf;
    for (int i = 0; i < n; ++i) {
      const double gumbel = -std::log(-std::log(unif_rand()));   // unif_rand() in (0,1)
      const double v = resid[i] * resid[i] * inv_2s + gumbel;
      if (v > best) { best = v; k = i; }
    }
    folds[iter] = k + 1;

    // ---- Step 2: theta | k  ~  p(theta | y_{-k})  via HMC -------------------
    double U_cur = potential_fold(Xp, yp, theta.data(), n, d, sigma2, tau2, k,
                                  resid.data(), g_cur.data());

    for (int h = 0; h < n_hmc; ++h) {
      for (int j = 0; j < d; ++j) {
        p[j]       = R::rnorm(0.0, 1.0);
        theta_p[j] = theta[j];
        p_p[j]     = p[j] - 0.5 * epsilon * g_cur[j];
      }

      double U_prop = R_NaReal;
      for (int l = 0; l < L; ++l) {
        for (int j = 0; j < d; ++j) theta_p[j] += epsilon * p_p[j];
        U_prop = potential_fold(Xp, yp, theta_p.data(), n, d, sigma2, tau2, k,
                                resid.data(), g.data());
        if (l != L - 1) {
          for (int j = 0; j < d; ++j) p_p[j] -= epsilon * g[j];
        }
      }
      for (int j = 0; j < d; ++j) p_p[j] -= 0.5 * epsilon * g[j];

      double K_cur = 0.0, K_prop = 0.0;
      for (int j = 0; j < d; ++j) {
        K_cur  += 0.5 * p[j] * p[j];
        K_prop += 0.5 * p_p[j] * p_p[j];
      }
      const double log_alpha = -(U_prop + K_prop) + (U_cur + K_cur);

      if (std::log(R::runif(0.0, 1.0)) < log_alpha) {   // NaN -> reject
        for (int j = 0; j < d; ++j) {
          theta[j] = theta_p[j];
          g_cur[j] = g[j];
        }
        U_cur = U_prop;
        acc += 1.0;
      }
    }

    for (int j = 0; j < d; ++j) samples(iter, j) = theta[j];

    if ((iter & 255) == 0) Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["samples"] = samples,
    _["acc"]     = acc / ((double) niters * n_hmc),   // HMC acceptance rate
    _["folds"]   = folds
  );
}
