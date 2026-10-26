// =============================================================================
// gibbs_loocv_ess.cpp  — component-wise elliptical slice within Gibbs
//
// Same target as gibbs_loocv.cpp, the LOOCV mixture proposal
//
//     q_mix(theta, phi, k)  ∝  p(theta | phi) p(phi) * p(y_{-k} | theta),   k = 1,...,n
//
// under the REGULARISED HORSESHOE prior (Piironen & Vehtari, 2017):
//
//   theta_0 ~ N(0, 100)                                   [intercept]
//   theta_j | lambda_j, tau, c ~ N(0, 1 / P_j),  P_j = 1/(tau^2 lambda_j^2) + 1/c^2
//   lambda_j ~ C+(0,1),  tau ~ C+(0, tau0),  c^2 ~ InvGamma(df/2, df*s^2/2)
//
// MODEL M1: Logistic link     p(y=1 | eta) = 1 / (1 + exp(-eta))
// MODEL M2: Cloglog link      p(y=1 | eta) = 1 - exp(-exp(eta))
//
// Gibbs iteration (HMC replaced by component-wise elliptical slice sampling):
//   1. Fold draw    k | theta            exact, Gumbel-max                [block]
//   2. Coefficients theta_j | theta_{-j}, phi, k   for j = 0,...,D
//      Univariate elliptical slice sampling (Murray, Adams & MacKay, 2010)
//      with the conditional prior N(0, 1/P_j) as the Gaussian factor and
//      prod_{i != k} p(y_i | eta_i) as the likelihood factor.  The linear
//      predictor eta is updated in O(n) per proposal.
//   3. Horseshoe    phi | theta
//      3a. lambda_1..lambda_D are conditionally independent given
//          (theta, tau, c^2), so the loop is one block update     [block]
//      3b. tau,  3c. joint (tau, lambda) shift,  3d. c^2   — slice sampling
//   4. Interweaving (ASIS; Yu & Meng, 2011) for tau:
//      z_j = theta_j * sqrt(P_j) held fixed, tau | z, lambda, c^2, theta_0, k
//      updated by slice sampling (it now sees the likelihood), then
//      theta_j = z_j / sqrt(P_j(tau)).
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace {

enum LinkType { LINK_LOGISTIC, LINK_CLOGLOG };

// ---------------------------------------------------------------------------
//  Log-likelihood of one observation
// ---------------------------------------------------------------------------
inline double softplus(double x) {
  return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

// log p(y | eta);  cloglog y=1 term is log(1 - exp(-exp(eta))) = log(-expm1(-u))
inline double loglik_obs(double y, double eta, LinkType link) {
  if (link == LINK_LOGISTIC) return y * eta - softplus(eta);
  const double u = std::exp(eta);
  return (y > 0.5) ? std::log(-std::expm1(-u)) : -u;
}

// ---------------------------------------------------------------------------
//  Slice sampler (Neal 2003: stepping-out + shrinkage)
// ---------------------------------------------------------------------------
template <class F>
inline double slice_sample(double x0, F logf, double w, int max_steps = 30) {
  const double logy = logf(x0) - exp_rand();
  double Lb = x0 - w * unif_rand(), Rb = Lb + w;
  int J = (int)std::floor(max_steps * unif_rand()), Kk = max_steps - 1 - J;
  while (J--  > 0 && logy < logf(Lb)) Lb -= w;
  while (Kk-- > 0 && logy < logf(Rb)) Rb += w;
  for (int it = 0; it < 500; ++it) {
    const double x1 = Lb + unif_rand() * (Rb - Lb);
    if (logy < logf(x1)) return x1;
    if (x1 < x0) Lb = x1; else Rb = x1;
  }
  return x0;
}

// ---------------------------------------------------------------------------
//  Horseshoe hyperparameter update  phi | theta  (identical to gibbs_loocv.cpp)
// ---------------------------------------------------------------------------
inline void update_horseshoe(
    arma::vec& ul, double& lt, double& lc2,
    const std::vector<double>& theta,
    double lt0, double half_sdf, double sdf_s2_half, double global_df,
    int D, arma::vec& a_j, arma::vec& b2_j,
    std::vector<double>& pprec, double tau2_int)
{
  const double ic2 = std::exp(-lc2);
  const double e2t = std::exp(-2.0 * lt);

  // 3a. Local scales lambda_j (conditionally independent: block update)
  for (int j = 0; j < D; ++j) {
    const double b2 = theta[1 + j] * theta[1 + j];
    auto logf = [&](double u) {
      if (u < -30.0 || u > 30.0) return R_NegInf;
      const double prec = std::exp(-2.0 * u) * e2t + ic2;
      return 0.5 * std::log(prec) - 0.5 * b2 * prec - softplus(2.0 * u) + u;
    };
    ul[j] = slice_sample(ul[j], logf, 1.5);
  }
  for (int j = 0; j < D; ++j) {
    a_j[j]  = std::exp(-2.0 * ul[j]);
    b2_j[j] = theta[1 + j] * theta[1 + j];
  }

  // 3b. Global scale tau
  {
    const double ic2c = std::exp(-lc2);
    auto logf = [&](double t) {
      if (t < -40.0 || t > 40.0) return R_NegInf;
      const double e2 = std::exp(-2.0 * t);
      double s = 0.0;
      for (int j = 0; j < D; ++j) {
        const double prec = a_j[j] * e2 + ic2c;
        s += 0.5 * std::log(prec) - 0.5 * b2_j[j] * prec;
      }
      return s - 0.5 * (global_df + 1.0) *
             std::log1p(std::exp(2.0 * (t - lt0)) / global_df) + t;
    };
    lt = slice_sample(lt, logf, 1.0);
  }

  // 3c. Joint scale move (tau * lambda_j preserved)
  {
    auto logf = [&](double sh) {
      const double tn = lt + sh;
      if (tn < -40.0 || tn > 40.0) return R_NegInf;
      double v = -0.5 * (global_df + 1.0) *
                 std::log1p(std::exp(2.0 * (tn - lt0)) / global_df) + tn;
      for (int j = 0; j < D; ++j) {
        const double un = ul[j] - sh;
        if (un < -30.0 || un > 30.0) return R_NegInf;
        v += -softplus(2.0 * un) + un;
      }
      return v;
    };
    const double sh = slice_sample(0.0, logf, 1.0);
    lt += sh;
    for (int j = 0; j < D; ++j) { ul[j] -= sh; a_j[j] = std::exp(-2.0 * ul[j]); }
  }

  // 3d. Slab variance c^2
  {
    const double e2t2 = std::exp(-2.0 * lt);
    auto logf = [&](double w) {
      if (w < -40.0 || w > 40.0) return R_NegInf;
      const double iw = std::exp(-w);
      double s = 0.0;
      for (int j = 0; j < D; ++j) {
        const double prec = a_j[j] * e2t2 + iw;
        s += 0.5 * std::log(prec) - 0.5 * b2_j[j] * prec;
      }
      return s - half_sdf * w - sdf_s2_half * iw;
    };
    lc2 = slice_sample(lc2, logf, 1.0);
  }

  // 3e. Refresh prior precision vector
  const double e2t_new = std::exp(-2.0 * lt), ic2_new = std::exp(-lc2);
  pprec[0] = 1.0 / tau2_int;
  for (int j = 0; j < D; ++j)
    pprec[1 + j] = a_j[j] * e2t_new + ic2_new;
}

// ---------------------------------------------------------------------------
//  Step 4: non-centred (interweaving) update of tau
//
//  With z_j = theta_j sqrt(P_j) and (lambda, c^2, theta_0, k) fixed,
//      theta_j(t) = z_j / sqrt(exp(-2 u_j - 2 t) + exp(-v)),   t = log tau,
//  and z ~ N(0, I) a priori does not depend on t, so
//      log q(t | z, .) = sum_{i != k} ll_i(eta_i(t)) + log C+(e^t; tau0) + t.
//  Each evaluation costs O(n D).
// ---------------------------------------------------------------------------
inline void asis_tau_update(
    std::vector<double>& theta, const arma::vec& ul, double& lt, double lc2,
    double lt0, double global_df, int D, int n,
    const double* Xp, const double* yp, int k, LinkType link,
    std::vector<double>& pprec, double tau2_int,
    std::vector<double>& z, std::vector<double>& eta0, std::vector<double>& eta_buf)
{
  for (int j = 0; j < D; ++j) z[j] = theta[1 + j] * std::sqrt(pprec[1 + j]);
  for (int i = 0; i < n; ++i) eta0[i] = Xp[i] * theta[0];
  const double ic2 = std::exp(-lc2);

  auto logf = [&](double t) {
    if (t < -40.0 || t > 40.0) return R_NegInf;
    const double e2t = std::exp(-2.0 * t);
    std::copy(eta0.begin(), eta0.end(), eta_buf.begin());
    for (int j = 0; j < D; ++j) {
      const double th = z[j] / std::sqrt(std::exp(-2.0 * ul[j]) * e2t + ic2);
      if (th == 0.0) continue;
      const double* col = Xp + (std::size_t)(1 + j) * n;
      for (int i = 0; i < n; ++i) eta_buf[i] += col[i] * th;
    }
    double s = 0.0;
    for (int i = 0; i < n; ++i) if (i != k) s += loglik_obs(yp[i], eta_buf[i], link);
    if (!std::isfinite(s)) return R_NegInf;
    return s - 0.5 * (global_df + 1.0) *
           std::log1p(std::exp(2.0 * (t - lt0)) / global_df) + t;
  };
  lt = slice_sample(lt, logf, 1.0);

  const double e2t = std::exp(-2.0 * lt);
  pprec[0] = 1.0 / tau2_int;
  for (int j = 0; j < D; ++j) {
    pprec[1 + j]  = std::exp(-2.0 * ul[j]) * e2t + ic2;
    theta[1 + j]  = z[j] / std::sqrt(pprec[1 + j]);
  }
}

// ---------------------------------------------------------------------------
//  Step 2: one univariate elliptical slice update of theta_j
//
//  Full conditional  ∝  N(theta_j; 0, 1/P_j) * exp{ sum_{i != k} ll_i(eta_i) }.
//  Proposal on the ellipse  theta_j' = theta_j cos(a) + nu sin(a),
//  nu ~ N(0, 1/P_j).  Only the change in log-likelihood is needed:
//      Delta(theta_j') = sum_{i != k} [ ll(eta_i + x_ij (theta_j' - theta_j)) - ll_i ].
//  On acceptance eta and the cached ll are updated in place.
//  Returns the number of likelihood evaluations used.
// ---------------------------------------------------------------------------
inline int ess_update_coord(
    int j, int k, int n,
    const double* col, const double* y,
    double* theta_j, double pprec_j,
    double* eta, double* ll, double* eta_new, double* ll_new,
    LinkType link)
{
  const double cur = *theta_j;
  const double nu  = norm_rand() / std::sqrt(pprec_j);
  const double log_u = std::log(unif_rand());

  double a     = 2.0 * M_PI * unif_rand();
  double a_min = a - 2.0 * M_PI, a_max = a;
  int n_eval = 0;

  for (int it = 0; it < 1000; ++it) {
    const double prop  = cur * std::cos(a) + nu * std::sin(a);
    const double delta = prop - cur;
    double dll = 0.0;
    for (int i = 0; i < n; ++i) {
      eta_new[i] = eta[i] + col[i] * delta;
      ll_new[i]  = loglik_obs(y[i], eta_new[i], link);
      if (i != k) dll += ll_new[i] - ll[i];
    }
    ++n_eval;
    if (std::isfinite(dll) && dll > log_u) {
      *theta_j = prop;
      std::copy(eta_new, eta_new + n, eta);
      std::copy(ll_new,  ll_new  + n, ll);
      return n_eval;
    }
    if (a < 0.0) a_min = a; else a_max = a;
    a = a_min + (a_max - a_min) * unif_rand();
  }
  return n_eval;  // bracket collapsed onto the current point: keep theta_j
}

// =============================================================================
//  Core ESS-within-Gibbs runner
// =============================================================================
List run_gibbs_loocv_rhs_ess(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& theta_init,
    double tau0,
    double slab_scale,
    double slab_df,
    int    niters,
    int    n_sweeps,
    double tau_init,
    double c2_init,
    bool   interweave_tau,
    LinkType link_type
) {
  const int n = X.n_rows, d = X.n_cols;
  const int D = d - 1;
  if ((int)y.n_elem != n || (int)theta_init.n_elem != d)
    stop("Dimension mismatch between X, y, and theta_init.");
  if (D < 1) stop("Need at least one penalised coefficient (d >= 2).");
  if (n_sweeps < 1) stop("n_sweeps must be >= 1.");

  const double* Xp = X.memptr();
  const double* yp = y.memptr();

  const double tau2_int    = 100.0;
  const double global_df   = 1.0;
  const double lt0         = std::log(tau0);
  const double half_sdf    = 0.5 * slab_df;
  const double sdf_s2_half = 0.5 * slab_df * slab_scale * slab_scale;

  // State
  std::vector<double> theta(d);
  for (int j = 0; j < d; ++j) theta[j] = theta_init[j];

  arma::vec ul(D, arma::fill::zeros), a_j(D), b2_j(D);
  double lt  = (tau_init > 0.0) ? std::log(tau_init) : lt0;
  double lc2 = (c2_init  > 0.0) ? std::log(c2_init)  : std::log(slab_scale * slab_scale);

  std::vector<double> pprec(d);
  {
    const double e2t = std::exp(-2.0 * lt), ic2 = std::exp(-lc2);
    pprec[0] = 1.0 / tau2_int;
    for (int j = 0; j < D; ++j) pprec[1 + j] = std::exp(-2.0 * ul[j]) * e2t + ic2;
  }

  // Linear predictor and per-observation log-likelihood caches
  std::vector<double> eta(n, 0.0), ll(n), eta_new(n), ll_new(n);
  std::vector<double> z_buf(D), eta0_buf(n);

  // Return structures
  NumericMatrix samples(niters, d);
  IntegerVector folds(niters);
  NumericVector tau_tr(niters), c2_tr(niters);
  double total_eval = 0.0;

  for (int it = 0; it < niters; ++it) {

    // Recompute eta and ll from scratch once per iteration (no drift)
    std::fill(eta.begin(), eta.end(), 0.0);
    for (int j = 0; j < d; ++j) {
      if (theta[j] == 0.0) continue;
      const double* col = Xp + (std::size_t)j * n;
      for (int i = 0; i < n; ++i) eta[i] += col[i] * theta[j];
    }
    for (int i = 0; i < n; ++i) ll[i] = loglik_obs(yp[i], eta[i], link_type);

    // -----------------------------------------------------------------------
    // Step 1: Fold draw k | theta  (Gumbel-max on -ll)
    // -----------------------------------------------------------------------
    int k = 0;
    double best = R_NegInf;
    for (int i = 0; i < n; ++i) {
      const double v = -ll[i] - std::log(-std::log(unif_rand()));
      if (v > best) { best = v; k = i; }
    }
    folds[it] = k + 1;

    // -----------------------------------------------------------------------
    // Step 2: theta_j | theta_{-j}, phi, k  by univariate ESS, j = 0..D
    // -----------------------------------------------------------------------
    for (int sw = 0; sw < n_sweeps; ++sw) {
      for (int j = 0; j < d; ++j) {
        total_eval += ess_update_coord(
          j, k, n, Xp + (std::size_t)j * n, yp,
          &theta[j], pprec[j],
          eta.data(), ll.data(), eta_new.data(), ll_new.data(),
          link_type);
      }
    }

    // -----------------------------------------------------------------------
    // Step 3: Horseshoe phi | theta
    // -----------------------------------------------------------------------
    update_horseshoe(ul, lt, lc2, theta,
                     lt0, half_sdf, sdf_s2_half, global_df,
                     D, a_j, b2_j, pprec, tau2_int);

    // -----------------------------------------------------------------------
    // Step 4: interweaving (non-centred) update of tau
    // -----------------------------------------------------------------------
    if (interweave_tau)
      asis_tau_update(theta, ul, lt, lc2, lt0, global_df, D, n,
                      Xp, yp, k, link_type, pprec, tau2_int,
                      z_buf, eta0_buf, eta_new);

    for (int j = 0; j < d; ++j) samples(it, j) = theta[j];
    tau_tr[it] = std::exp(lt);
    c2_tr[it]  = std::exp(lc2);

    if ((it & 127) == 0) Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["samples"]        = samples,
    _["folds"]          = folds,
    _["tau"]            = tau_tr,
    _["c2"]             = c2_tr,
    _["evals_per_coord"] = total_eval / ((double)niters * n_sweeps * d)
  );
}

} // anonymous namespace

// =============================================================================
//  Exported R function
// =============================================================================

// [[Rcpp::export]]
List gibbs_loocv_rhs_ess_cpp(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& theta_init,
    double tau0,
    double slab_scale = 2.0,
    double slab_df    = 4.0,
    int    niters     = 10000,
    std::string link  = "logit",
    int    n_sweeps   = 1,
    double tau_init   = -1.0,
    double c2_init    = -1.0,
    bool   interweave_tau = true
) {
  LinkType lt;
  if (link == "logit" || link == "logistic") lt = LINK_LOGISTIC;
  else if (link == "cloglog")                lt = LINK_CLOGLOG;
  else stop("Unknown link: '" + link + "'. Use 'logit' or 'cloglog'.");

  return run_gibbs_loocv_rhs_ess(
    X, y, theta_init, tau0, slab_scale, slab_df, niters,
    n_sweeps, tau_init, c2_init, interweave_tau, lt
  );
}
