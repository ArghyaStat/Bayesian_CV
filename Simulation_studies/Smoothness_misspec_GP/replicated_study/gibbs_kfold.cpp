// [[Rcpp::depends(RcppArmadillo)]]
// -----------------------------------------------------------------------------
// Univariate Vecchia GP regression with nugget: Gibbs sampler on the K-fold mixture
//
//   y = X beta + w + eps,   w ~ GP(0, sigma2 K_{phi,nu}),   eps ~ N(0, tau2 I),
//   s2 = sigma2 + tau2,     r = sigma2 / s2,
//   Cov(y | theta) = s2 C_{phi,r},   C_{phi,r} = r K_{phi,nu} + (1 - r) I,
//   beta | s2 ~ N(0, s2 V),   s2 ~ IG(a, b),   phi ~ U(0, b_phi),   r ~ U(0, 1).
//
//   q_mix(theta, U = k)  ∝  p(theta) p(y | theta) / p(y_{D_k} | y_{-D_k}, theta)
//                        =  p(theta) p(y_{-D_k} | theta),   theta = (beta, s2, phi, r).
//
// One sweep:
//   U | theta                     P(U = k | theta) ∝ 1 / p(y_{D_k} | y_{-D_k}, theta)
//   z | theta, U = k              z ~ p(y_{D_k} | y_{-D_k}, theta)        (augmentation)
//   phi | r, y~, then r | phi, y~ random-walk Metropolis-Hastings on p(phi, r | y~),
//                                 (beta, s2) integrated out, in the unconstrained
//                                 (u, v) = (logit(phi / b_phi), logit(r))
//   (s2, beta) | phi, r, y~       s2 ~ IG (beta integrated out), beta | s2 ~ N
// where y~ = (y_{-D_k}, z). Given U = k the last three steps form a blocked Gibbs
// cycle on p(theta, z | y_{-D_k}) and hence leave p(theta | y_{-D_k}) invariant.
//
// All densities are taken under the Vecchia joint N(X beta, s2 (U'U)^{-1}) of the
// response, where row i of U e is (e_i - sum_j B_ij e_{N(i)_j}) / sqrt(d_i) and
// (B_i, d_i) are the kriging weights and variance of C_{phi,r}.
// -----------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>

using namespace Rcpp;

// ------------------------------------------------------------------ structures

// Vecchia conditioning sets (0-based, location i excluded) with indices into a
// table of unique pairwise distances, so that each (phi, r) requires one Matern
// evaluation per distinct distance only.
struct Vecchia {
  int N;
  std::vector<std::vector<int> > nb;
  std::vector<std::vector<int> > d22;  // lower triangle of K(N(i), N(i)), row-wise
  std::vector<std::vector<int> > d21;  // K(N(i), i)
  std::vector<double> ud;
};

// Covariance parameters (phi, r) with the Vecchia whitener U and sum_i log d_i
struct CovState {
  double phi;
  double r;
  arma::sp_mat U;
  double sum_log_d;
  bool valid;
};

// Held-out block D_k: selection matrix P_k (N x n_k) and Cholesky factor of
// Q_kk = (U P_k)'(U P_k), the unscaled precision of y_{D_k} | y_{-D_k}
struct Fold {
  arma::uvec idx;
  arma::sp_mat P;
};

struct FoldFactor {
  arma::sp_mat UH;
  arma::mat R;
  double logdet;
};

// Conjugate quantities given (phi, r, y~):
//   P = (UX)'(UX) + V^{-1} = R'R,  c = R^{-T} (UX)'(Uy),  ssr = ||Uy||^2 - ||c||^2,
//   log_marg = log p(y~ | phi, r) + const
//            = -sum_i log d_i / 2 - log|R| - (a + N/2) log(b + ssr/2)
struct Collapsed {
  arma::mat R;
  arma::vec c;
  double ssr;
  double log_marg;
};

// Triangular solves with an upper Cholesky factor R: R' x = b and R x = b
static arma::vec solve_Rt(const arma::mat& R, const arma::vec& b) {
  return arma::solve(arma::trimatl(R.t()), b, arma::solve_opts::fast);
}

static arma::vec solve_R(const arma::mat& R, const arma::vec& b) {
  return arma::solve(arma::trimatu(R), b, arma::solve_opts::fast);
}

// ------------------------------------------------------------- Vecchia helpers

// Matern correlation, fields::Matern(dist, range = phi, smoothness = nu):
//   rho(h) = 2^{1-nu} / Gamma(nu) (h / phi)^nu K_nu(h / phi)
// [[Rcpp::export]]
double matern(double dist, double phi, double nu) {
  if (dist == 0.0) return 1.0;
  const double x = dist / phi;

  // Half-integer nu = p + 1/2:
  //   rho = exp(-x) p!/(2p)! sum_{i=0}^p (p+i)! / (i! (p-i)!) (2x)^{p-i}
  const double p = nu - 0.5;
  if (p >= 0.0 && p == std::floor(p)) {
    const int ip = (int)p;
    double s = 0.0;
    for (int i = 0; i <= ip; ++i)
      s += std::exp(std::lgamma(ip + i + 1.0) - std::lgamma(i + 1.0) -
                    std::lgamma(ip - i + 1.0) + (ip - i) * std::log(2.0 * x));
    return std::exp(-x + std::lgamma(ip + 1.0) - std::lgamma(2.0 * ip + 1.0)) * s;
  }

  return std::exp((1.0 - nu) * std::log(2.0) - std::lgamma(nu) + nu * std::log(x) +
                  std::log(R::bessel_k(x, nu, 1.0)));
}

static Vecchia build_structure(const IntegerMatrix& NN, const NumericMatrix& locs) {
  const int N = NN.nrow(), mp1 = NN.ncol(), d = locs.ncol();
  Vecchia V;
  V.N = N;
  V.nb.assign(N, std::vector<int>());
  V.d22.assign(N, std::vector<int>());
  V.d21.assign(N, std::vector<int>());

  std::map<long long, int> key;
  auto dist_index = [&](int a, int b) {
    double ss = 0.0;
    for (int c = 0; c < d; ++c) {
      const double df = locs(a, c) - locs(b, c);
      ss += df * df;
    }
    const double dist = std::sqrt(ss);
    const long long k = std::llround(dist * 1e10);
    auto it = key.find(k);
    if (it != key.end()) return it->second;
    const int id = (int)V.ud.size();
    V.ud.push_back(dist);
    key[k] = id;
    return id;
  };

  for (int i = 0; i < N; ++i) {
    bool self = false;
    for (int c = 0; c < mp1; ++c) {
      const int v = NN(i, c);
      if (v == NA_INTEGER) continue;
      if (v - 1 == i) { self = true; continue; }
      if (v - 1 > i) stop("NNarray row %d contains a later location.", i + 1);
      V.nb[i].push_back(v - 1);
    }
    if (!self) stop("NNarray row %d does not contain the location itself.", i + 1);
    const std::vector<int>& nbi = V.nb[i];
    for (std::size_t a = 0; a < nbi.size(); ++a) {
      for (std::size_t b = 0; b <= a; ++b) V.d22[i].push_back(dist_index(nbi[a], nbi[b]));
      V.d21[i].push_back(dist_index(nbi[a], i));
    }
  }
  return V;
}

// Vecchia factor of C_{phi,r}: correlation r rho(h) for h > 0 and 1 at h = 0
static CovState cov_state(const Vecchia& V, double phi, double r, double nu) {
  const int N = V.N;
  arma::vec kv(V.ud.size());
  for (std::size_t u = 0; u < V.ud.size(); ++u)
    kv(u) = V.ud[u] == 0.0 ? 1.0 : r * matern(V.ud[u], phi, nu);

  arma::uword nnz = N;
  for (int i = 0; i < N; ++i) nnz += V.nb[i].size();
  arma::umat pos(2, nnz);
  arma::vec val(nnz);
  arma::uword t = 0;
  double sum_log_d = 0.0;

  for (int i = 0; i < N; ++i) {
    const int n1 = V.nb[i].size();
    double d = 1.0;
    arma::vec x;
    if (n1 > 0) {
      arma::mat K22(n1, n1);
      arma::vec K21(n1);
      int c = 0;
      for (int a = 0; a < n1; ++a) {
        for (int b = 0; b <= a; ++b) K22(a, b) = K22(b, a) = kv(V.d22[i][c++]);
        K21(a) = kv(V.d21[i][a]);
      }
      arma::mat R;
      if (!arma::chol(R, K22, "upper")) return CovState{phi, r, arma::sp_mat(), 0.0, false};
      x = solve_R(R, solve_Rt(R, K21));
      d = 1.0 - arma::dot(x, K21);
    }
    if (!(d > 0.0) || !std::isfinite(d)) return CovState{phi, r, arma::sp_mat(), 0.0, false};
    const double s = 1.0 / std::sqrt(d);
    pos(0, t) = i; pos(1, t) = i; val(t++) = s;
    for (int a = 0; a < n1; ++a) {
      pos(0, t) = i; pos(1, t) = V.nb[i][a]; val(t++) = -x(a) * s;
    }
    sum_log_d += std::log(d);
  }
  return CovState{phi, r, arma::sp_mat(pos, val, N, N), sum_log_d, true};
}

static std::vector<Fold> build_folds(const IntegerVector& fold, int K) {
  const int N = fold.size();
  std::vector<std::vector<arma::uword> > members(K);
  for (int i = 0; i < N; ++i) {
    if (fold[i] < 1 || fold[i] > K) stop("fold must take values in 1..K.");
    members[fold[i] - 1].push_back(i);
  }
  std::vector<Fold> F(K);
  for (int k = 0; k < K; ++k) {
    const arma::uword n = members[k].size();
    if (n == 0) stop("Fold %d is empty.", k + 1);
    F[k].idx = arma::uvec(members[k]);
    arma::umat pos(2, n);
    for (arma::uword h = 0; h < n; ++h) { pos(0, h) = F[k].idx(h); pos(1, h) = h; }
    F[k].P = arma::sp_mat(pos, arma::ones<arma::vec>(n), N, n);
  }
  return F;
}

static std::vector<FoldFactor> fold_factors(const CovState& cs, const std::vector<Fold>& F) {
  std::vector<FoldFactor> out(F.size());
  for (std::size_t k = 0; k < F.size(); ++k) {
    out[k].UH = cs.U * F[k].P;
    out[k].R = arma::chol(arma::mat(out[k].UH.t() * out[k].UH), "upper");
    out[k].logdet = 2.0 * arma::sum(arma::log(out[k].R.diag()));
  }
  return out;
}

// y_{D_k} | y_{-D_k} ~ N(X_k beta + m_k, s2 Q_kk^{-1}),
//   m_k = -Q_kk^{-1} Q_{k,-k} e_{-k},   Q_kk^{-1} Q_{k,-k} e_{-k} = R^{-1} a_k,
//   a_k = R^{-T} (U P_k)'(U e - U P_k e_k)
static arma::vec fold_shift(const arma::vec& Ue, const arma::vec& eH, const FoldFactor& f) {
  const arma::vec rhs = f.UH.t() * (Ue - f.UH * eH);
  return solve_Rt(f.R, rhs);
}

// log p(y_{D_l} | y_{-D_l}, theta), l = 1..K, at the observed y
static arma::vec fold_lik(const arma::vec& e, double s2, const CovState& cs,
                          const std::vector<Fold>& F, const std::vector<FoldFactor>& FF) {
  const arma::vec Ue = cs.U * e;
  arma::vec out(F.size());
  for (std::size_t l = 0; l < F.size(); ++l) {
    const arma::vec eH = e.elem(F[l].idx);
    const arma::vec z = FF[l].R * eH + fold_shift(Ue, eH, FF[l]);
    const double nH = eH.n_elem;
    out(l) = -0.5 * (nH * std::log(2.0 * M_PI * s2) - FF[l].logdet + arma::dot(z, z) / s2);
  }
  return out;
}

// ---------------------------------------------------------------- Gibbs steps

// P(U = k | theta) ∝ exp(-fold_ll_k), Gumbel-max draw
static int update_U(const arma::vec& fold_ll) {
  int k = 0;
  double best = R_NegInf;
  for (arma::uword l = 0; l < fold_ll.n_elem; ++l) {
    const double v = -fold_ll(l) - std::log(-std::log(unif_rand()));
    if (v > best) { best = v; k = l; }
  }
  return k;
}

// z ~ p(y_{D_k} | y_{-D_k}, theta); returns y~ with the D_k block replaced
static arma::vec impute_heldout(const arma::vec& y, const arma::vec& mu, double s2,
                                const CovState& cs, const Fold& Fk, const FoldFactor& f) {
  const arma::vec e = y - mu;
  const arma::vec eH = e.elem(Fk.idx);
  const arma::vec a = fold_shift(cs.U * e, eH, f);
  arma::vec xi(eH.n_elem);
  for (arma::uword h = 0; h < xi.n_elem; ++h) xi(h) = norm_rand();
  arma::vec out = y;
  out.elem(Fk.idx) = mu.elem(Fk.idx) + solve_R(f.R, std::sqrt(s2) * xi - a);
  return out;
}

static Collapsed collapse(const arma::vec& y_aug, const arma::mat& X, const CovState& cs,
                          const arma::mat& Vinv, double a_sigma, double b_sigma) {
  const arma::mat UX = cs.U * X;
  const arma::vec Uy = cs.U * y_aug;
  Collapsed cl;
  cl.R = arma::chol(UX.t() * UX + Vinv, "upper");
  cl.c = solve_Rt(cl.R, UX.t() * Uy);
  cl.ssr = std::max(arma::dot(Uy, Uy) - arma::dot(cl.c, cl.c), 0.0);
  cl.log_marg = -0.5 * cs.sum_log_d - arma::sum(arma::log(cl.R.diag())) -
    (a_sigma + 0.5 * y_aug.n_elem) * std::log(b_sigma + 0.5 * cl.ssr);
  return cl;
}

// ------------------------------------------- Metropolis-Hastings for (phi, r)

static inline double softplus(double x) { return std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x))); }
static inline double expit(double x) { return 1.0 / (1.0 + std::exp(-x)); }
static inline double logit(double p) { return std::log(p) - std::log1p(-p); }

// Target of the (phi, r) update for the current y~
struct Target {
  const arma::vec& y_aug;
  const arma::mat& X;
  const Vecchia& V;
  const arma::mat& Vinv;
  double a_sigma, b_sigma, nu, b_phi;
};

// Point (u, v) with its Vecchia factor, conjugate quantities and log target
//   lp = log p(phi, r | y~) + log|d(phi, r) / d(u, v)| + const,
//   log d expit(x) / dx = -softplus(-x) - softplus(x)
struct Point {
  double u, v;
  CovState cs;
  Collapsed cl;
  double lp;
};

static Point eval_point(double u, double v, const Target& tg) {
  Point pt{u, v, cov_state(tg.V, tg.b_phi * expit(u), expit(v), tg.nu), Collapsed(), R_NegInf};
  if (!pt.cs.valid || !(pt.cs.phi > 0.0) || !(pt.cs.r > 0.0) || !(pt.cs.r < 1.0)) return pt;
  pt.cl = collapse(tg.y_aug, tg.X, pt.cs, tg.Vinv, tg.a_sigma, tg.b_sigma);
  pt.lp = pt.cl.log_marg - softplus(-u) - softplus(u) - softplus(-v) - softplus(v);
  return pt;
}

// Gaussian random walk on u (j = 0) or v (j = 1) with variance tuning; the proposal
// is symmetric on the logit scale, so the acceptance ratio is exp(lp* - lp)
static bool update_logit(Point& cur, int j, double tuning, const Target& tg) {
  const double step = std::sqrt(tuning) * norm_rand();
  Point can = j == 0 ? eval_point(cur.u + step, cur.v, tg) : eval_point(cur.u, cur.v + step, tg);
  if (std::log(unif_rand()) < can.lp - cur.lp) {
    cur = std::move(can);
    return true;
  }
  return false;
}

// (s2, beta) | phi, r, y~:
//   s2 ~ IG(a + N/2, b + ssr/2),   beta | s2 ~ N(R^{-1} c, s2 (R'R)^{-1})
static void update_beta_s2(arma::vec& beta, double& s2, const Collapsed& cl, int N,
                           double a_sigma, double b_sigma) {
  s2 = 1.0 / R::rgamma(a_sigma + 0.5 * N, 1.0 / (b_sigma + 0.5 * cl.ssr));
  arma::vec xi(beta.n_elem);
  for (arma::uword j = 0; j < xi.n_elem; ++j) xi(j) = norm_rand();
  beta = solve_R(cl.R, cl.c + std::sqrt(s2) * xi);
}

// ---------------------------------------------------------------- main sampler

// [[Rcpp::export]]
List gibbs_kfold(const arma::vec& y, const arma::mat& X, const NumericMatrix& locs,
                 const IntegerMatrix& NNarray, const IntegerVector& fold, int K,
                 double nu, double phi_init, double r_init, const arma::vec& beta_init,
                 double s2_init, const arma::mat& V_prior, double a_sigma, double b_sigma,
                 double b_phi, double tuning_phi, double tuning_r, int niters,
                 int fixed_fold = 0) {

  const int N = y.n_elem, p = X.n_cols;
  if ((int)X.n_rows != N || locs.nrow() != N || NNarray.nrow() != N || fold.size() != N ||
      (int)beta_init.n_elem != p)
    stop("Dimension mismatch between inputs.");
  if (phi_init <= 0.0 || phi_init >= b_phi) stop("phi_init must lie in (0, b_phi).");
  if (r_init <= 0.0 || r_init >= 1.0) stop("r_init must lie in (0, 1).");
  if (fixed_fold < 0 || fixed_fold > K) stop("fixed_fold must lie in 0..K.");

  const Vecchia V = build_structure(NNarray, locs);
  const std::vector<Fold> F = build_folds(fold, K);
  const arma::mat Vinv = arma::inv_sympd(V_prior);

  // initial state theta^(0)
  arma::vec beta = beta_init;
  double s2 = s2_init;
  CovState cs = cov_state(V, phi_init, r_init, nu);
  if (!cs.valid) stop("Initial (phi, r) gives a singular Vecchia factor.");
  std::vector<FoldFactor> FF = fold_factors(cs, F);
  arma::vec fold_ll = fold_lik(y - X * beta, s2, cs, F, FF);

  // storage: phi, r, beta, sigma2 = r s2, tau2 = (1 - r) s2
  NumericMatrix fold_ll_out(niters, K);
  NumericMatrix theta_out(niters, p + 4);
  IntegerVector k_out(niters);
  CharacterVector cn(p + 4);
  cn[0] = "phi";
  cn[1] = "r";
  for (int j = 0; j < p; ++j) cn[j + 2] = "beta" + std::to_string(j);
  cn[p + 2] = "sigma2";
  cn[p + 3] = "tau2";
  colnames(theta_out) = cn;

  int acc_phi = 0, acc_r = 0;
  const int step_10pct = std::max(1, niters / 10);
  const auto t_start = std::chrono::steady_clock::now();

  for (int iter = 0; iter < niters; ++iter) {

    // U | theta
    const int k = fixed_fold > 0 ? fixed_fold - 1 : update_U(fold_ll);

    // theta | U = k, invariant for p(theta | y_{-D_k})
    const arma::vec y_aug = impute_heldout(y, X * beta, s2, cs, F[k], FF[k]);

    const Target tg{y_aug, X, V, Vinv, a_sigma, b_sigma, nu, b_phi};
    Point pt = eval_point(logit(cs.phi / b_phi), logit(cs.r), tg);
    const bool a_phi = update_logit(pt, 0, tuning_phi, tg);
    const bool a_r   = update_logit(pt, 1, tuning_r, tg);
    acc_phi += a_phi;
    acc_r   += a_r;
    if (a_phi || a_r) {
      cs = pt.cs;
      FF = fold_factors(cs, F);
    }

    update_beta_s2(beta, s2, pt.cl, N, a_sigma, b_sigma);

    // log p(y_{D_l} | y_{-D_l}, theta), l = 1..K, at the new theta
    fold_ll = fold_lik(y - X * beta, s2, cs, F, FF);

    // record
    for (int l = 0; l < K; ++l) fold_ll_out(iter, l) = fold_ll(l);
    k_out[iter] = k + 1;
    theta_out(iter, 0) = cs.phi;
    theta_out(iter, 1) = cs.r;
    for (int j = 0; j < p; ++j) theta_out(iter, j + 2) = beta(j);
    theta_out(iter, p + 2) = cs.r * s2;
    theta_out(iter, p + 3) = (1.0 - cs.r) * s2;

    if ((iter & 255) == 0) checkUserInterrupt();
    if ((iter + 1) % step_10pct == 0 || (iter + 1) == niters) {
      const double mins = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count() / 60.0;
      Rcout << "Progress: " << (int)((iter + 1) * 100.0 / niters) << "% ("
            << (iter + 1) << " / " << niters << ") - Elapsed: "
            << std::setprecision(2) << std::fixed << mins << " mins\n";
    }
  }

  const double rate_phi = niters > 0 ? (double)acc_phi / niters : NA_REAL;
  const double rate_r   = niters > 0 ? (double)acc_r / niters : NA_REAL;

  return List::create(
    _["samples"]     = theta_out,
    _["theta"]       = theta_out,
    _["fold_loglik"] = fold_ll_out,
    _["acc"]         = rate_phi,
    _["acc_phi"]     = rate_phi,
    _["acc_r"]       = rate_r,
    _["folds"]       = k_out,
    _["state"]       = List::create(
      _["beta"]   = beta,
      _["sigma2"] = cs.r * s2,
      _["tau2"]   = (1.0 - cs.r) * s2,
      _["phi"]    = cs.phi,
      _["r"]      = cs.r
    )
  );
}
