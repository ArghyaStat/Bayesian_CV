// [[Rcpp::depends(RcppArmadillo)]]
// Model M2 (separate covariances): Gibbs sampler on the mixture q_mix(theta, U) for K-fold CV
// U | theta (fold indicator), then phi, sig2, B, W | U = k invariant for p(theta | y_{-D_k})
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

// Vecchia coefficients B and conditional variances r for a given phi
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

// Sufficient quantities for the conjugate (B, sig2) | W, phi update
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

// log N_Vecchia(W | M, K_phi, diag(sig2)):  columns are independent GPs
static double vecchia_ll_sep(const arma::mat& W, const arma::mat& M, const Vecchia& V,
                             const Coeff& C, const arma::vec& sig2) {
  const int N = W.n_rows, q = W.n_cols;
  const arma::mat E = W - M;
  double ll = 0.0;
  for (int j = 0; j < q; ++j) {
    const double s2 = sig2(j);
    for (int i = 0; i < N; ++i) {
      const std::vector<int>& nbi = V.nb[i];
      const std::vector<double>& Bi = C.B[i];
      double res = E(i, j);
      int b = 0;
      for (std::size_t t = 0; t < nbi.size(); ++t) {
        if ((int)t == V.self[i]) continue;
        res -= Bi[b++] * E(nbi[t], j);
      }
      const double ri = C.r[i];
      ll += -0.5 * (std::log(2.0 * M_PI * ri * s2) + res * res / (ri * s2));
    }
  }
  return ll;
}

// Vecchia-whitened products U W and U X
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

// P(U = k | theta) ∝ exp(-log p(y_{D_k} | W_{D_k}))   (Gumbel-max draw)
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


// phi | W, B, sig2  (uniform prior on (0, b_phi), truncated-normal RW proposal)
static PhiState update_phi(const PhiState& cur, const arma::mat& M, const arma::mat& W,
                           const arma::vec& sig2, const NumericMatrix& locs,
                           const Vecchia& V, double nu, double b_phi, double tuning) {
  const double s = std::sqrt(tuning);
  
  double can;
  do { can = cur.value + s * norm_rand(); } while (can <= 0.0 || can >= b_phi);
  
  const Coeff  Ccan     = vecchia_coeff_fill(locs, V, can, nu);
  const double post_cur = vecchia_ll_sep(W, M, V, cur.C, sig2);
  const double post_can = vecchia_ll_sep(W, M, V, Ccan,  sig2);
  
  // Hastings correction for the truncated proposal
  const double Zcur = R::pnorm((b_phi - cur.value) / s, 0, 1, 1, 0) - R::pnorm(-cur.value / s, 0, 1, 1, 0);
  const double Zcan = R::pnorm((b_phi - can)       / s, 0, 1, 1, 0) - R::pnorm(-can       / s, 0, 1, 1, 0);
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


// sig2_j | W_j, phi ~ IG(a + N/2, b + c_j/2), beta_j integrated out; Binomial scale fixed
static arma::vec update_sig2(const RegStats& rs, const arma::mat& Mprior,
                             const arma::mat& cholP, double a_prior, double b_prior,
                             const std::vector<int>& fam, const std::vector<double>& fixed_sd) {
  const int N = rs.UW.n_rows, q = rs.UW.n_cols;
  arma::vec sig2(q);
  for (int j = 0; j < q; ++j) {
    if (fam[j] == 2) { sig2(j) = fixed_sd[j] * fixed_sd[j]; continue; }
    const arma::vec a = cholP * Mprior.col(j);
    const arma::vec b = rs.Rtilde * rs.Mpart.col(j);
    const double c = arma::dot(rs.UW.col(j), rs.UW.col(j)) + arma::dot(a, a) - arma::dot(b, b);
    sig2(j) = 1.0 / R::rgamma(a_prior + 0.5 * N, 1.0 / (b_prior + 0.5 * c));
  }
  return sig2;
}

// beta_j | sig2_j, W_j, phi ~ N(V_tilde Mpart_j, sig2_j V_tilde)
static arma::mat update_beta_sep(const RegStats& rs, const arma::vec& sig2) {
  const int p = rs.Mpart.n_rows, q = rs.Mpart.n_cols;
  const arma::mat Z = rnorm_mat(p, q);
  arma::mat beta(p, q);
  for (int j = 0; j < q; ++j)
    beta.col(j) = rs.Rtilde.t() * (rs.Rtilde * rs.Mpart.col(j) + std::sqrt(sig2(j)) * Z.col(j));
  return beta;
}

// Sparse Vecchia whitening matrix U, with precision Q = U'U
static arma::sp_mat vecchia_whitener(const Vecchia& V, const Coeff& C) {
  arma::uword nnz = V.N;
  for (const auto& b : C.B) nnz += b.size();
  arma::umat positions(2, nnz);
  arma::vec values(nnz);
  arma::uword t = 0;
  for (int i = 0; i < V.N; ++i) {
    if (!(C.r[i] > 0.0) || !std::isfinite(C.r[i]))
      stop("Nonpositive or nonfinite Vecchia conditional variance.");
    const double scale = 1.0 / std::sqrt(C.r[i]);
    positions(0, t) = i; positions(1, t) = i; values(t++) = scale;
    int b = 0;
    for (std::size_t a = 0; a < V.nb[i].size(); ++a) {
      if ((int)a == V.self[i]) continue;
      positions(0, t) = i; positions(1, t) = V.nb[i][a];
      values(t++) = -C.B[i][b++] * scale;
    }
  }
  return arma::sp_mat(positions, values, V.N, V.N);
}

// Held-out latent rows W_{D_k} | W_{-D_k}, B, sig2, phi
static arma::mat draw_heldout_sep(const arma::mat& W, const arma::mat& M,
                                  const arma::vec& sig2, const Vecchia& V,
                                  const Coeff& C, const arma::uvec& held) {
  if (held.is_empty()) return W;
  const arma::sp_mat U = vecchia_whitener(V, C); // Q = U' U
  // Extract held-out columns of U without forming a dense N x N matrix.
  arma::uword nnz = 0;
  for (arma::uword h = 0; h < held.n_elem; ++h)
    for (auto it = U.begin_col(held(h)); it != U.end_col(held(h)); ++it) ++nnz;
  arma::umat positions(2, nnz);
  arma::vec values(nnz);
  arma::uword t = 0;
  for (arma::uword h = 0; h < held.n_elem; ++h) {
    for (auto it = U.begin_col(held(h)); it != U.end_col(held(h)); ++it) {
      positions(0, t) = it.row(); positions(1, t) = h; values(t++) = *it;
    }
  }
  const arma::sp_mat UH(positions, values, V.N, held.n_elem);
  const arma::mat QHH(UH.t() * UH);
  const arma::mat R = arma::chol(QHH, "upper");
  arma::mat E = W - M;
  E.rows(held).zeros();
  const arma::mat rhs = UH.t() * (U * E); // Q_HT (W_T - M_T)
  const arma::mat mean = M.rows(held) - arma::solve(arma::trimatu(R),
                                arma::solve(arma::trimatl(R.t()), rhs));
  const arma::mat noise = arma::solve(arma::trimatu(R),
                                      rnorm_mat(held.n_elem, W.n_cols));
  arma::mat out = W;
  arma::mat scaled_noise = noise;
  scaled_noise.each_row() %= arma::sqrt(sig2).t();
  out.rows(held) = mean + scaled_noise;
  return out;
}

// W | B, sig2, phi, y_{-D_k}: elliptical slice on training sites, then held-out draw
static arma::mat update_W_sep(const arma::mat& W_cur, const arma::mat& M,
                              const arma::vec& sig2, const Vecchia& V,
                              const Coeff& C, const arma::mat& Y,
                              const std::vector<int>& fam, const std::vector<double>& par,
                              const std::vector<char>& train) {
  const int N = W_cur.n_rows, q = W_cur.n_cols;
  std::vector<arma::uword> tv, hv;
  for (int i = 0; i < N; ++i) (train[i] ? tv : hv).push_back(i);
  const arma::uvec training(tv), held(hv);
  arma::mat W = W_cur;
  if (!training.is_empty()) {
    // Prior draw restricted to training sites gives the slice directions
    const arma::mat directions = vecchia_prior_draw(rnorm_mat(N, q), V, C);
    for (int j = 0; j < q; ++j) {
      arma::vec mu = M.col(j);
      const double Scond = sig2(j);
      if (!(Scond > 0.0) || !std::isfinite(Scond))
        stop("Nonpositive or nonfinite conditional response variance.");
      const double sc = std::sqrt(Scond);
      auto ll_col = [&](const arma::vec& col) {
        double ll = 0.0;
        for (arma::uword t = 0; t < training.n_elem; ++t) {
          const arma::uword i = training(t);
          ll += loglik1(fam[j], par[j], Y(i, j), col(i));
        }
        return ll;
      };
      const arma::vec cur = W.col(j);
      const double ll_cur = ll_col(cur);
      if (!std::isfinite(ll_cur)) stop("Nonfinite training log likelihood in latent update.");
      const double logy = ll_cur + std::log(unif_rand());
      double angle = unif_rand() * 2.0 * M_PI;
      double lower = angle - 2.0 * M_PI, upper = angle;
      arma::vec candidate = cur;
      // Slice-bracket search; continues until acceptance
      unsigned int attempts = 0;
      while (true) {
        if ((attempts++ & 63u) == 0u) Rcpp::checkUserInterrupt();
        for (arma::uword t = 0; t < training.n_elem; ++t) {
          const arma::uword i = training(t);
          candidate(i) = mu(i) + (cur(i) - mu(i)) * std::cos(angle)
            + sc * directions(i, j) * std::sin(angle);
        }
        const double ll = ll_col(candidate);
        if (std::isnan(ll)) stop("NaN training log likelihood in slice search.");
        if (ll > logy) {
          for (arma::uword t = 0; t < training.n_elem; ++t)
            W(training(t), j) = candidate(training(t));
          break;
        }
        if (angle < 0.0) lower = angle; else upper = angle;
        angle = lower + (upper - lower) * unif_rand();
      }
    }
  }
  // Refresh held-out latent rows from their Gaussian conditional
  return draw_heldout_sep(W, M, sig2, V, C, held);
}

// ---------------------------------------------------------------- main sampler

// [[Rcpp::export]]
List gibbs_sep(
    const NumericMatrix& Y, const NumericMatrix& X, const NumericMatrix& locs,
    const IntegerMatrix& NNarray, const IntegerVector& fold, int K,
    const CharacterVector& family, const NumericVector& family_par,
    double nu, double phi_init, const NumericMatrix& beta_init,
    const NumericVector& sig2_init, const NumericMatrix& W_init,
    const NumericMatrix& M_prior, const NumericMatrix& V_prior,
    double a_prior, double b_prior,
    double b_phi, double tuning_phi,
    int niters) {
  
  const int N = Y.nrow(), q = Y.ncol(), p = X.ncol();
  if (X.nrow() != N || locs.nrow() != N || NNarray.nrow() != N || fold.size() != N ||
      W_init.nrow() != N || W_init.ncol() != q || sig2_init.size() != q) {
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
  const arma::mat precP  = arma::inv_sympd(as<arma::mat>(V_prior));
  const arma::mat cholP  = arma::chol(precP, "upper");
  const Vecchia   V      = build_structure(NNarray);
  
  // initial state theta^(0)
  arma::mat W    = as<arma::mat>(W_init);
  arma::mat beta = as<arma::mat>(beta_init);
  arma::vec sig2 = as<arma::vec>(sig2_init);
  PhiState  phi{phi_init, vecchia_coeff_fill(locs, V, phi_init, nu), false};   // phi.value, phi.C
  arma::vec fold_ll = fold_lik(Ya, W, fold0, K, fam, par);
  
  std::vector<double> fixed_sd(q);
  for (int j = 0; j < q; ++j) fixed_sd[j] = std::sqrt(sig2(j));
  
  // storage
  const int ntheta = 1 + p * q + q;
  NumericMatrix fold_ll_out(niters, K);
  NumericMatrix theta_out(niters, ntheta);
  IntegerVector k_out(niters);
  
  int acc_phi = 0;
  const int step_10pct = std::max(1, niters / 10);
  const auto t_start = std::chrono::steady_clock::now();
  
  for (int iter = 0; iter < niters; ++iter) {
    
    // U | theta
    const int k = update_U(fold_ll);
    const std::vector<char> train = train_mask(fold0, k);
    
    // theta | U = k, invariant for p(theta | y_{-D_k})
    phi      = update_phi(phi, Xa * beta, W, sig2, locs, V, nu, b_phi, tuning_phi);
    acc_phi += phi.accepted;
    
    const RegStats rs = reg_stats(W, Xa, V, phi.C, precP, Mprior);
    sig2     = update_sig2(rs, Mprior, cholP, a_prior, b_prior, fam, fixed_sd);
    beta     = update_beta_sep(rs, sig2);
    
    W        = update_W_sep(W, Xa * beta, sig2, V, phi.C, Ya, fam, par, train);
    
    // Fold log-likelihoods log p(y_{D_l} | W_{D_l}) at the new theta
    fold_ll  = fold_lik(Ya, W, fold0, K, fam, par);
    
    // record
    for (int l = 0; l < K; ++l) fold_ll_out(iter, l) = fold_ll(l);
    k_out[iter] = k + 1;
    int c = 0;
    theta_out(iter, c++) = phi.value;
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < p; ++i) theta_out(iter, c++) = beta(i, j);
    for (int j = 0; j < q; ++j) theta_out(iter, c++) = sig2(j);
    
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
      _["W"]    = W,
      _["beta"] = beta,
      _["sig2"] = sig2,
      _["phi"]  = phi.value
    )
  );
}
