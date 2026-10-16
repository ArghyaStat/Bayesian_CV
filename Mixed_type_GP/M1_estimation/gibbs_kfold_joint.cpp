// [[Rcpp::depends(RcppArmadillo)]]
// -----------------------------------------------------------------------------
// Model M1 (joint covariance): Gibbs sampler on the augmented mixture
//
//   q_mix(theta, U = k)  ∝  p(theta) p(y | theta) p(y_{D_k} | W_{D_k})^{-1},
//   theta = (B, Sigma, phi, W_{-D_k}, W_{D_k}),
//
// alternating the two full conditionals of eljpd_derivation.pdf:
//
//   (Eq. 7)  P(U = k | theta)  ∝  p(y_{D_k} | W_{D_k})^{-1},          k = 1..K
//   (Eq. 8)  q_mix(theta | U = k)  ∝  p(theta | y_{-D_k})
//              ∝ p(y_{-D_k} | W_{-D_k}) p(W | B, Sigma, phi) p(B | Sigma) p(Sigma) p(phi)
//              = p(W_{D_k} | W_{-D_k}, B, Sigma, phi) p(W_{-D_k}, B, Sigma, phi | y_{-D_k}),
//            i.e. the training likelihood enters exactly once.
//
// Given U = k, theta is refreshed by the blocked Gibbs cycle (Algorithm 1, steps 3-6)
//   phi   | W, B, Sigma                  Metropolis-Hastings (Vecchia GP density)
//   Sigma | W, phi                       inverse-Wishart, B integrated out
//   B     | Sigma, W, phi                matrix normal
//   W     | B, Sigma, phi, y_{-D_k}      elliptical slice sampling, training likelihood only
// each of which leaves p(theta | y_{-D_k}) invariant.
//
// W_{D_k} is part of the chain state theta.  The fold log-likelihoods
// log p(y_{D_l} | W_{D_l}) used in Eq. 7 and returned for the SNIS estimator are
// evaluated at the CURRENT W, never at a freshly simulated latent field.
//
// Every update below is a pure function: it takes the current state by const
// reference and RETURNS the new value, e.g. W = update_W(...),
// phi = update_phi(...), Sigma = update_Sigma(...), beta = update_beta(...),
// k = update_U(...).  Nothing is modified in place or stored globally.
// -----------------------------------------------------------------------------
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <iomanip>

using namespace Rcpp;

// ---------------------------------------------------------------- structures

// Vecchia neighbour structure (0-based); nb[i] contains i itself at position self[i]
struct Vecchia {
  int N;
  std::vector<std::vector<int> > nb;
  std::vector<int> self;
};

// Vecchia regression coefficients and conditional variances for a given phi:
//   W_i | W_{N(i)} ~ N( sum_a B[i][a] W_{N(i)_a},  r[i] * Sigma )
struct Coeff {
  std::vector<std::vector<double> > B;
  std::vector<double> r;
};

// State of the range parameter: phi together with its Vecchia coefficients
struct PhiState {
  double value;
  Coeff  C;
  bool   accepted;
};

// Sufficient quantities for the conjugate (B, Sigma) | W, phi update
struct RegStats {
  arma::mat UW;      // U W   (Vecchia-whitened latent field), N x q
  arma::mat Mpart;   // X'U'U W + V^{-1} M_prior,             p x q
  arma::mat Rtilde;  // upper chol of V_tilde = (X'U'U X + V^{-1})^{-1}
};

// ------------------------------------------------------------ Vecchia helpers

static Vecchia build_structure(const IntegerMatrix& NN) {
  const int N = NN.nrow(), mp1 = NN.ncol();
  Vecchia V;
  V.N = N;
  V.nb.assign(N, std::vector<int>());
  V.self.assign(N, -1);
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < mp1; ++c) {
      const int v = NN(i, c);
      if (v == NA_INTEGER) continue;
      V.nb[i].push_back(v - 1);
    }
    for (std::size_t t = 0; t < V.nb[i].size(); ++t) {
      if (V.nb[i][t] == i) V.self[i] = (int)t;
    }
    if (V.self[i] < 0) stop("NNarray row %d does not contain the location itself.", i + 1);
  }
  return V;
}

// [[Rcpp::export]]
double matern(double dist, double phi, double nu) {
  if (dist == 0.0) return 1.0;
  const double kappa = std::sqrt(2.0 * nu) / phi;
  const double part1 = std::pow(2.0, 1.0 - nu) / std::tgamma(nu);
  const double part2 = kappa * dist;
  const double out = std::log(part1) + nu * std::log(part2) +
                     std::log(R::bessel_k(part2, nu, 1.0));
  return std::exp(out);
}

static Coeff vecchia_coeff_fill(const NumericMatrix& locs, const Vecchia& V,
                           double phi, double nu) {
  const int N = V.N, d = locs.ncol();
  Coeff C;
  C.B.resize(N);
  C.r.resize(N);

  for (int i = 0; i < N; ++i) {
    const std::vector<int>& nbi = V.nb[i];
    const int n = (int)nbi.size();
    if (n == 1) {
      C.B[i].clear();
      C.r[i] = 1.0;
      continue;
    }

    const int n1 = n - 1, s = V.self[i];
    std::vector<int> oth;
    oth.reserve(n1);
    for (int t = 0; t < n; ++t) if (t != s) oth.push_back(nbi[t]);

    arma::mat K22(n1, n1);
    arma::vec K21(n1);
    for (int a = 0; a < n1; ++a) {
      for (int b = 0; b <= a; ++b) {
        double ss = 0.0;
        for (int c = 0; c < d; ++c) {
          const double df = locs(oth[a], c) - locs(oth[b], c);
          ss += df * df;
        }
        const double val = matern(std::sqrt(ss), phi, nu);
        K22(a, b) = val;
        K22(b, a) = val;
      }
      double ss = 0.0;
      for (int c = 0; c < d; ++c) {
        const double df = locs(oth[a], c) - locs(i, c);
        ss += df * df;
      }
      K21(a) = matern(std::sqrt(ss), phi, nu);
    }

    const arma::mat R = arma::chol(K22, "upper");
    const arma::vec x = arma::solve(arma::trimatu(R),
                                    arma::solve(arma::trimatl(R.t()), K21));
    C.B[i] = arma::conv_to<std::vector<double> >::from(x);
    C.r[i] = 1.0 - arma::dot(x, K21);
  }
  return C;
}

// log N_Vecchia(W | M, K_phi, Sigma)   (matrix-normal, rows via Vecchia, cols via Sigma)
static double vecchia_ll(const arma::mat& W, const arma::mat& M, const Vecchia& V,
                         const Coeff& C, const arma::mat& cholSigma) {
  const int N = W.n_rows, q = W.n_cols;
  const arma::mat E = W - M;
  double logdetS = 0.0;
  for (int j = 0; j < q; ++j) logdetS += 2.0 * std::log(cholSigma(j, j));

  std::vector<double> res(q), z(q);
  double ll = 0.0;
  for (int i = 0; i < N; ++i) {
    const std::vector<int>& nbi = V.nb[i];
    const std::vector<double>& Bi = C.B[i];
    for (int j = 0; j < q; ++j) res[j] = E(i, j);
    int b = 0;
    for (std::size_t t = 0; t < nbi.size(); ++t) {
      if ((int)t == V.self[i]) continue;
      const double w = Bi[b++];
      for (int j = 0; j < q; ++j) res[j] -= w * E(nbi[t], j);
    }
    // z = L^{-1} res  with Sigma = L L', L = cholSigma'
    for (int j = 0; j < q; ++j) {
      double t = res[j];
      for (int l = 0; l < j; ++l) t -= cholSigma(l, j) * z[l];
      z[j] = t / cholSigma(j, j);
    }
    double quad = 0.0;
    for (int j = 0; j < q; ++j) quad += z[j] * z[j];
    const double ri = C.r[i];
    ll += -0.5 * (q * std::log(2.0 * M_PI * ri) + logdetS + quad / ri);
  }
  return ll;
}

// Vecchia-whitened products U W and U X, where row i of U A is
//   (A_i - sum_a B[i][a] A_{N(i)_a}) / sqrt(r_i)
struct VecchiaProducts {
  arma::mat UW;
  arma::mat UX;
};

static VecchiaProducts vecchia_products(const arma::mat& W, const arma::mat& X,
                                        const Vecchia& V, const Coeff& C) {
  const int N = W.n_rows, q = W.n_cols, p = X.n_cols;
  VecchiaProducts P;
  P.UW.set_size(N, q);
  P.UX.set_size(N, p);
  for (int i = 0; i < N; ++i) {
    const std::vector<int>& nbi = V.nb[i];
    const std::vector<double>& Bi = C.B[i];
    const double is = 1.0 / std::sqrt(C.r[i]);
    for (int j = 0; j < q; ++j) P.UW(i, j) = W(i, j);
    for (int j = 0; j < p; ++j) P.UX(i, j) = X(i, j);
    int b = 0;
    for (std::size_t t = 0; t < nbi.size(); ++t) {
      if ((int)t == V.self[i]) continue;
      const double w = Bi[b++];
      for (int j = 0; j < q; ++j) P.UW(i, j) -= w * W(nbi[t], j);
      for (int j = 0; j < p; ++j) P.UX(i, j) -= w * X(nbi[t], j);
    }
    for (int j = 0; j < q; ++j) P.UW(i, j) *= is;
    for (int j = 0; j < p; ++j) P.UX(i, j) *= is;
  }
  return P;
}

// Draw from N_Vecchia(0, K_phi) independently for each column of Z ~ N(0, I)
static arma::mat vecchia_prior_draw(const arma::mat& Z, const Vecchia& V, const Coeff& C) {
  const int N = Z.n_rows, q = Z.n_cols;
  arma::mat out(N, q);
  for (int i = 0; i < N; ++i) {          // neighbours precede i in the ordering
    const std::vector<int>& nbi = V.nb[i];
    const std::vector<double>& Bi = C.B[i];
    const double sr = std::sqrt(C.r[i]);
    for (int j = 0; j < q; ++j) out(i, j) = sr * Z(i, j);
    int b = 0;
    for (std::size_t t = 0; t < nbi.size(); ++t) {
      if ((int)t == V.self[i]) continue;
      const double w = Bi[b++];
      for (int j = 0; j < q; ++j) out(i, j) += w * out(nbi[t], j);
    }
  }
  return out;
}

static arma::mat rnorm_mat(int n, int m) {
  arma::mat Z(n, m);
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i) Z(i, j) = norm_rand();
  return Z;
}

// -------------------------------------------------------- data-level likelihood

static inline double loglik1(int fam, double par, double y, double eta) {
  switch (fam) {
  case 0: return R::dnorm(y, eta, std::sqrt(par), 1);
  case 1: return R::dpois(y, std::exp(eta), 1);
  case 2: {
    const double pr = (eta >= 0.0) ? 1.0 / (1.0 + std::exp(-eta))
                                   : std::exp(eta) / (1.0 + std::exp(eta));
    return R::dbinom(y, par, pr, 1);
  }
  case 3: {
    const double mu = std::exp(eta);
    return R::dgamma(y, par, mu / par, 1);
  }
  case 4: return R::dnbinom_mu(y, par, std::exp(eta), 1);
  }
  return R_NegInf;
}

// log p(y_{D_l} | W_{D_l}) for l = 1..K, evaluated at the CURRENT latent field W
static arma::vec fold_lik(const arma::mat& Y, const arma::mat& W,
                             const std::vector<int>& fold0, int K,
                             const std::vector<int>& fam, const std::vector<double>& par) {
  const int N = Y.n_rows, q = Y.n_cols;
  arma::vec out(K, arma::fill::zeros);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int j = 0; j < q; ++j) s += loglik1(fam[j], par[j], Y(i, j), W(i, j));
    out(fold0[i]) += s;
  }
  return out;
}

static std::vector<char> train_mask(const std::vector<int>& fold0, int k) {
  std::vector<char> train(fold0.size());
  for (std::size_t i = 0; i < fold0.size(); ++i) train[i] = (fold0[i] != k);
  return train;
}

// ------------------------------------------------------------ Gibbs updates

// Eq. 7:  P(U = k | theta) ∝ exp(-log p(y_{D_k} | W_{D_k}))   (Gumbel-max draw)
static int update_U(const arma::vec& fold_ll) {
  const int K = fold_ll.n_elem;
  int k = 0;
  double best = R_NegInf;
  for (int l = 0; l < K; ++l) {
    const double v = -fold_ll(l) - std::log(-std::log(unif_rand()));
    if (v > best) { best = v; k = l; }
  }
  return k;
}

// phi | W, B, Sigma  (uniform prior on (0, b_phi), truncated-normal RW proposal)
static PhiState update_phi(const PhiState& cur, const arma::mat& M, const arma::mat& W,
                           const arma::mat& cholSigma, const NumericMatrix& locs,
                           const Vecchia& V, double nu, double b_phi, double tuning) {
  const double s = std::sqrt(tuning);

  double can;
  do { can = cur.value + s * norm_rand(); } while (can <= 0.0 || can >= b_phi);

  const Coeff  Ccan     = vecchia_coeff_fill(locs, V, can, nu);
  const double post_cur = vecchia_ll(W, M, V, cur.C, cholSigma);
  const double post_can = vecchia_ll(W, M, V, Ccan,  cholSigma);

  // Hastings correction for the truncated proposal
  const double Zcur = R::pnorm((b_phi - cur.value) / s, 0, 1, 1, 0) - R::pnorm(-cur.value / s, 0, 1, 1, 0);
  const double Zcan = R::pnorm((b_phi - can)     / s, 0, 1, 1, 0) - R::pnorm(-can     / s, 0, 1, 1, 0);
  const double log_r = (post_can - post_cur) + std::log(Zcur) - std::log(Zcan);

  if (std::log(unif_rand()) < log_r) return PhiState{can, Ccan, true};
  return PhiState{cur.value, cur.C, false};
}

static RegStats reg_stats(const arma::mat& W, const arma::mat& X,
                                        const Vecchia& V, const Coeff& C,
                                        const arma::mat& precP, const arma::mat& Mprior) {
  RegStats rs;
  const VecchiaProducts P = vecchia_products(W, X, V, C);
  rs.UW = P.UW;
  const arma::mat& UX = P.UX;
  const arma::mat V_tilde = arma::inv_sympd(UX.t() * UX + precP);
  rs.Rtilde = arma::chol(V_tilde, "upper");
  rs.Mpart  = UX.t() * rs.UW + precP * Mprior;
  return rs;
}

static arma::mat riwish_cpp(double df, const arma::mat& S) {
  const int q = S.n_rows;
  const arma::mat R = arma::chol(arma::inv_sympd(S), "upper");
  arma::mat A(q, q, arma::fill::zeros);
  for (int i = 0; i < q; ++i) {
    A(i, i) = std::sqrt(R::rchisq(df - i));
    for (int j = i + 1; j < q; ++j) A(i, j) = norm_rand();
  }
  const arma::mat AR = A * R;
  return arma::inv_sympd(AR.t() * AR);
}

// Sigma | W, phi  (B integrated out); Binomial margins keep their fixed scale
static arma::mat update_Sigma(const RegStats& rs, const arma::mat& Mprior,
                              const arma::mat& cholP, const arma::mat& Sprior,
                              double df_prior, const std::vector<int>& fam,
                              const std::vector<double>& fixed_sd) {
  const int N = rs.UW.n_rows, q = rs.UW.n_cols;
  const arma::mat a = cholP * Mprior;
  const arma::mat b = rs.Rtilde * rs.Mpart;
  arma::mat S_tilde = Sprior + rs.UW.t() * rs.UW + a.t() * a - b.t() * b;
  S_tilde = 0.5 * (S_tilde + S_tilde.t());

  const arma::mat Sstar = riwish_cpp(df_prior + N, S_tilde);
  const arma::vec dstar = arma::sqrt(Sstar.diag());
  arma::mat Rstar = Sstar.each_col() / dstar;
  Rstar.each_row() /= dstar.t();

  arma::vec dnew = dstar;
  for (int j = 0; j < q; ++j) if (fam[j] == 2) dnew(j) = fixed_sd[j];

  arma::mat Sigma = Rstar;
  Sigma.each_col() %= dnew;
  Sigma.each_row() %= dnew.t();
  return Sigma;
}

// B | Sigma, W, phi  ~  MN(V_tilde Mpart, V_tilde, Sigma)
static arma::mat update_beta(const RegStats& rs, const arma::mat& cholSigma) {
  const arma::mat Z = rnorm_mat(rs.Mpart.n_rows, rs.Mpart.n_cols);
  return rs.Rtilde.t() * (rs.Rtilde * rs.Mpart + Z * cholSigma);
}

// Eq. 8, latent block:  W | B, Sigma, phi, y_{-D_k}
//   ∝ p(y_{-D_k} | W_{-D_k}) p(W_{D_k} | W_{-D_k}, .) p(W_{-D_k} | .)
// Elliptical slice sampling column by column (W_j | W_{-j}); the likelihood
// only includes training sites, so W_{D_k} moves under its GP conditional.
static arma::mat update_W(const arma::mat& W_cur, const arma::mat& M, const arma::mat& Sigma,
                          const Vecchia& V, const Coeff& C, const arma::mat& Y,
                          const std::vector<int>& fam, const std::vector<double>& par,
                          const std::vector<char>& train) {
  const int N = W_cur.n_rows, q = W_cur.n_cols;
  arma::mat W = W_cur;                                   // local copy, returned
  const arma::mat nu = vecchia_prior_draw(rnorm_mat(N, q), V, C);

  arma::vec mu(N), cand(N);

  for (int j = 0; j < q; ++j) {
    // column j | other columns:  N(mu, Scond * K_phi)
    double Scond = Sigma(j, j);
    arma::vec cvec;
    arma::uvec minus;
    if (q > 1) {
      minus.set_size(q - 1);
      int c = 0;
      for (int l = 0; l < q; ++l) if (l != j) minus(c++) = l;
      const arma::vec Smj = Sigma.submat(minus, arma::uvec({(arma::uword)j}));
      cvec  = arma::solve(Sigma.submat(minus, minus), Smj);
      Scond = Sigma(j, j) - arma::dot(Smj, cvec);
    }
    for (int i = 0; i < N; ++i) {
      double m = M(i, j);
      for (int c = 0; c < q - 1; ++c) m += (W(i, minus(c)) - M(i, minus(c))) * cvec(c);
      mu(i) = m;
    }
    const double sc = std::sqrt(Scond);

    auto ll_col = [&](const arma::vec& col) {
      double s = 0.0;
      for (int i = 0; i < N; ++i)
        if (train[i]) s += loglik1(fam[j], par[j], Y(i, j), col(i));
      return s;
    };

    const arma::vec cur = W.col(j);
    const double logy = ll_col(cur) + std::log(unif_rand());

    double gam  = unif_rand() * 2.0 * M_PI;
    double gmin = gam - 2.0 * M_PI, gmax = gam;
    for (int it = 0; it < 1000; ++it) {
      const double cg = std::cos(gam), sg = std::sin(gam);
      for (int i = 0; i < N; ++i)
        cand(i) = mu(i) + (cur(i) - mu(i)) * cg + sc * nu(i, j) * sg;
      if (ll_col(cand) > logy) { W.col(j) = cand; break; }
      if (gam < 0.0) gmin = gam; else gmax = gam;
      gam = gmin + (gmax - gmin) * unif_rand();
    }
  }
  return W;
}

// ---------------------------------------------------------------- main sampler

// [[Rcpp::export]]
List gibbs_joint(
    const NumericMatrix& Y, const NumericMatrix& X, const NumericMatrix& locs,
    const IntegerMatrix& NNarray, const IntegerVector& fold, int K,
    const CharacterVector& family, const NumericVector& family_par,
    double nu, double phi_init, const NumericMatrix& beta_init,
    const NumericMatrix& Sigma_init, const NumericMatrix& W_init,
    const NumericMatrix& M_prior, const NumericMatrix& V_prior,
    const NumericMatrix& S_prior, double df_prior,
    double b_phi, double tuning_phi,
    int niters) {

  const int N = Y.nrow(), q = Y.ncol(), p = X.ncol();
  if (X.nrow() != N || locs.nrow() != N || NNarray.nrow() != N || fold.size() != N ||
      W_init.nrow() != N || W_init.ncol() != q) {
    stop("Dimension mismatch between inputs.");
  }
  if (family.size() != q) stop("length(family) must equal ncol(Y).");
  if (phi_init <= 0.0 || phi_init >= b_phi) stop("phi_init must lie in (0, b_phi).");

  std::vector<int> fam(q);
  std::vector<double> par(q, 1.0);
  for (int j = 0; j < q; ++j) {
    const std::string f = as<std::string>(family[j]);
    if      (f == "Gaussian")          fam[j] = 0;
    else if (f == "Poisson")           fam[j] = 1;
    else if (f == "Binomial")          fam[j] = 2;
    else if (f == "Gamma")             fam[j] = 3;
    else if (f == "Negative-Binomial") fam[j] = 4;
    else stop("Unsupported family: " + f);
    if (family_par.size() == q && !NumericVector::is_na(family_par[j])) par[j] = family_par[j];
  }

  std::vector<int> fold0(N);
  for (int i = 0; i < N; ++i) {
    fold0[i] = fold[i] - 1;
    if (fold0[i] < 0 || fold0[i] >= K) stop("fold must take values in 1..K.");
  }

  // fixed data / prior quantities
  const arma::mat Ya     = as<arma::mat>(Y);
  const arma::mat Xa     = as<arma::mat>(X);
  const arma::mat Mprior = as<arma::mat>(M_prior);
  const arma::mat Sprior = as<arma::mat>(S_prior);
  const arma::mat precP  = arma::inv_sympd(as<arma::mat>(V_prior));
  const arma::mat cholP  = arma::chol(precP, "upper");
  const Vecchia   V      = build_structure(NNarray);

  // initial state theta^(0)
  arma::mat W         = as<arma::mat>(W_init);
  arma::mat beta      = as<arma::mat>(beta_init);
  arma::mat Sigma     = as<arma::mat>(Sigma_init);
  arma::mat cholSigma = arma::chol(Sigma, "upper");
  PhiState  phi{phi_init, vecchia_coeff_fill(locs, V, phi_init, nu), false};   // phi.value, phi.C
  arma::vec fold_ll   = fold_lik(Ya, W, fold0, K, fam, par);     // Step 1

  std::vector<double> fixed_sd(q);
  for (int j = 0; j < q; ++j) fixed_sd[j] = std::sqrt(Sigma(j, j));

  // storage
  const int ntheta = 1 + p * q + q * (q + 1) / 2;
  NumericMatrix fold_ll_out(niters, K);
  NumericMatrix theta_out(niters, ntheta);
  IntegerVector k_out(niters);

  int acc_phi = 0;
  const int step_10pct = std::max(1, niters / 10);
  const auto t_start = std::chrono::steady_clock::now();

  for (int iter = 0; iter < niters; ++iter) {

    // Step 2 (Eq. 7): U | theta
    const int k = update_U(fold_ll);
    const std::vector<char> train = train_mask(fold0, k);

    // Steps 3-6 (Eq. 8): theta | U = k, invariant for p(theta | y_{-D_k})
    phi       = update_phi(phi, Xa * beta, W, cholSigma, locs, V, nu, b_phi, tuning_phi);
    acc_phi  += phi.accepted;

    const RegStats rs = reg_stats(W, Xa, V, phi.C, precP, Mprior);
    Sigma     = update_Sigma(rs, Mprior, cholP, Sprior, df_prior, fam, fixed_sd);
    cholSigma = arma::chol(Sigma, "upper");
    beta      = update_beta(rs, cholSigma);

    W         = update_W(W, Xa * beta, Sigma, V, phi.C, Ya, fam, par, train);

    // Step 7: log p(y_{D_l} | W_{D_l}), l = 1..K, at the new theta
    fold_ll   = fold_lik(Ya, W, fold0, K, fam, par);

    // record
    for (int l = 0; l < K; ++l) fold_ll_out(iter, l) = fold_ll(l);
    k_out[iter] = k + 1;
    int c = 0;
    theta_out(iter, c++) = phi.value;
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < p; ++i) theta_out(iter, c++) = beta(i, j);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i <= j; ++i) theta_out(iter, c++) = Sigma(i, j);

    if ((iter & 15) == 0) Rcpp::checkUserInterrupt();
    if ((iter + 1) % step_10pct == 0 || (iter + 1) == niters) {
      const double mins = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count() / 60.0;
      Rcout << "Progress: " << (int)((iter + 1) * 100.0 / niters) << "% ("
            << (iter + 1) << " / " << niters << ") - Elapsed: "
            << std::setprecision(2) << std::fixed << mins << " mins\n";
    }
  }

  const double rate_phi = niters > 0 ? (double)acc_phi / niters : NA_REAL;

  return List::create(
    _["samples"]     = theta_out,
    _["theta"]       = theta_out,
    _["fold_loglik"] = fold_ll_out,
    _["acc"]         = rate_phi,
    _["acc_phi"]     = rate_phi,
    _["folds"]       = k_out,
    _["state"]       = List::create(
      _["W"]     = W,
      _["beta"]  = beta,
      _["Sigma"] = Sigma,
      _["phi"]   = phi.value
    )
  );
}
