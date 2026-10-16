// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>

using namespace Rcpp;

struct Vecchia {
  int N;
  std::vector<std::vector<int> > nb;
  std::vector<int> self;
};

struct Coeff {
  std::vector<std::vector<double> > B;
  std::vector<double> r;
};

static void build_structure(const IntegerMatrix& NN, Vecchia& V) {
  const int N = NN.nrow(), mp1 = NN.ncol();
  V.N = N;
  V.nb.assign(N, std::vector<int>());
  V.self.assign(N, -1);
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < mp1; ++c) {
      int v = NN(i, c);
      if (v == NA_INTEGER) continue;
      V.nb[i].push_back(v - 1);
    }
    for (std::size_t t = 0; t < V.nb[i].size(); ++t) {
      if (V.nb[i][t] == i) V.self[i] = (int)t;
    }
    if (V.self[i] < 0) stop("NNarray row %d does not contain the location itself.", i + 1);
  }
}

// [[Rcpp::export]]
double matern(double dist, double phi, double nu) {
  if (dist == 0.0) return 1.0;
  double kappa = std::sqrt(2.0 * nu) / phi;
  double part1 = std::pow(2.0, 1.0 - nu) / std::tgamma(nu);
  double part2 = kappa * dist;
  double out = std::log(part1) + nu * std::log(part2) + std::log(R::bessel_k(part2, nu, 1.0));
  return std::exp(out);
}

static void vecchia_coeff_fill(const NumericMatrix& locs, const Vecchia& V,
                               double phi, double nu, Coeff& C) {
  const int N = V.N, d = locs.ncol();
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
    for (int t = 0; t < n; ++t) {
      if (t != s) oth.push_back(nbi[t]);
    }

    arma::mat K22(n1, n1);
    arma::vec K21(n1);

    for (int a = 0; a < n1; ++a) {
      for (int b = 0; b <= a; ++b) {
        double ss = 0.0;
        for (int c = 0; c < d; ++c) {
          double df = locs(oth[a], c) - locs(oth[b], c);
          ss += df * df;
        }
        double val = matern(std::sqrt(ss), phi, nu);
        K22(a, b) = val;
        K22(b, a) = val;
      }
      double ss = 0.0;
      for (int c = 0; c < d; ++c) {
        double df = locs(oth[a], c) - locs(i, c);
        ss += df * df;
      }
      K21(a) = matern(std::sqrt(ss), phi, nu);
    }

    arma::mat R = arma::chol(K22, "upper");
    arma::vec y = arma::solve(arma::trimatl(R.t()), K21);
    arma::vec x = arma::solve(arma::trimatu(R), y);
    double dot = arma::dot(x, K21);

    std::vector<double> B_i(n1);
    for (int a = 0; a < n1; ++a) B_i[a] = x(a);

    C.B[i] = B_i;
    C.r[i] = 1.0 - dot;
  }
}

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

static void vecchia_products(const arma::mat& W, const arma::mat& X, const Vecchia& V,
                             const Coeff& C, arma::mat& UW, arma::mat& UX) {
  const int N = W.n_rows, q = W.n_cols, p = X.n_cols;
  UW.set_size(N, q);
  UX.set_size(N, p);
  for (int i = 0; i < N; ++i) {
    const std::vector<int>& nbi = V.nb[i];
    const std::vector<double>& Bi = C.B[i];
    const double is = 1.0 / std::sqrt(C.r[i]);
    for (int j = 0; j < q; ++j) UW(i, j) = W(i, j);
    for (int j = 0; j < p; ++j) UX(i, j) = X(i, j);
    int b = 0;
    for (std::size_t t = 0; t < nbi.size(); ++t) {
      if ((int)t == V.self[i]) continue;
      const double w = Bi[b++];
      for (int j = 0; j < q; ++j) UW(i, j) -= w * W(nbi[t], j);
      for (int j = 0; j < p; ++j) UX(i, j) -= w * X(nbi[t], j);
    }
    for (int j = 0; j < q; ++j) UW(i, j) *= is;
    for (int j = 0; j < p; ++j) UX(i, j) *= is;
  }
}

static void vecchia_prior_draw(const arma::mat& Z, const Vecchia& V, const Coeff& C,
                               arma::mat& nu) {
  const int N = Z.n_rows, q = Z.n_cols;
  nu.set_size(N, q);
  for (int i = 0; i < N; ++i) {
    const std::vector<int>& nbi = V.nb[i];
    const std::vector<double>& Bi = C.B[i];
    const double sr = std::sqrt(C.r[i]);
    for (int j = 0; j < q; ++j) nu(i, j) = sr * Z(i, j);
    int b = 0;
    for (std::size_t t = 0; t < nbi.size(); ++t) {
      if ((int)t == V.self[i]) continue;
      const double w = Bi[b++];
      for (int j = 0; j < q; ++j) nu(i, j) += w * nu(nbi[t], j);
    }
  }
}

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

// Out-of-sample Vecchia predictive log-likelihood for all folds (Model M2)
static void compute_all_fold_predictive_loglik_sep(
    const arma::mat& Y, const arma::mat& X, const arma::mat& W,
    const arma::mat& beta, const arma::vec& sig2,
    const NumericMatrix& locs, const List& test_nn_list,
    double phi, double nu, const std::vector<int>& fam, const std::vector<double>& par,
    std::vector<double>& out) {

  const int K = test_nn_list.size(), q = Y.n_cols, d = locs.ncol();
  const arma::mat M = X * beta;
  out.assign(K, 0.0);

  for (int k = 0; k < K; ++k) {
    List f_info = test_nn_list[k];
    IntegerVector t_idx = f_info["test_idx"];
    IntegerMatrix nn_idx = f_info["nn_idx"];
    const int n_test = t_idx.size(), m = nn_idx.ncol();

    arma::mat K22(m, m);
    arma::vec K21(m);
    double ll_fold = 0.0;

    for (int i = 0; i < n_test; ++i) {
      const int site_i = t_idx[i];

      for (int a = 0; a < m; ++a) {
        int n_a = nn_idx(i, a);
        for (int b = 0; b <= a; ++b) {
          int n_b = nn_idx(i, b);
          double ss = 0.0;
          for (int c = 0; c < d; ++c) {
            double df = locs(n_a, c) - locs(n_b, c);
            ss += df * df;
          }
          double val = matern(std::sqrt(ss), phi, nu);
          K22(a, b) = val; K22(b, a) = val;
        }
        double ss = 0.0;
        for (int c = 0; c < d; ++c) {
          double df = locs(n_a, c) - locs(site_i, c);
          ss += df * df;
        }
        K21(a) = matern(std::sqrt(ss), phi, nu);
      }

      arma::mat R = arma::chol(K22, "upper");
      arma::vec y = arma::solve(arma::trimatl(R.t()), K21);
      arma::vec x = arma::solve(arma::trimatu(R), y);
      double dot = arma::dot(x, K21);
      double r_i = 1.0 - dot;
      if (r_i < 1e-8) r_i = 1e-8;

      arma::rowvec mean_i = M.row(site_i);
      for (int a = 0; a < m; ++a) {
        int n_a = nn_idx(i, a);
        mean_i += x(a) * (W.row(n_a) - M.row(n_a));
      }

      for (int j = 0; j < q; ++j) {
        double w_draw = mean_i(j) + std::sqrt(r_i * sig2(j)) * norm_rand();
        ll_fold += loglik1(fam[j], par[j], Y(site_i, j), w_draw);
      }
    }
    out[k] = ll_fold;
  }
}

static void update_W_sep(arma::mat& W, const arma::mat& M, const arma::vec& sig2,
                         const Vecchia& V, const Coeff& C, const arma::mat& Y,
                         const std::vector<int>& fam, const std::vector<double>& par,
                         const std::vector<char>& train) {
  const int N = W.n_rows, q = W.n_cols;

  arma::mat Z(N, q);
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < N; ++i) Z(i, j) = norm_rand();
  }
  arma::mat nu;
  vecchia_prior_draw(Z, V, C, nu);

  arma::vec mu(N), cand(N);

  for (int j = 0; j < q; ++j) {
    const double sc = std::sqrt(sig2(j));
    for (int i = 0; i < N; ++i) mu(i) = M(i, j);

    auto ll_col = [&](const arma::vec& col) {
      double s = 0.0;
      for (int i = 0; i < N; ++i) {
        if (train[i]) s += loglik1(fam[j], par[j], Y(i, j), col(i));
      }
      return s;
    };

    arma::vec cur = W.col(j);
    const double ll_cur = ll_col(cur);
    const double logy = ll_cur + std::log(unif_rand());

    double gam = unif_rand() * 2.0 * M_PI;
    double gmin = gam - 2.0 * M_PI, gmax = gam;

    for (int it = 0; it < 1000; ++it) {
      const double cg = std::cos(gam), sg = std::sin(gam);
      for (int i = 0; i < N; ++i) {
        cand(i) = mu(i) + (cur(i) - mu(i)) * cg + sc * nu(i, j) * sg;
      }
      if (ll_col(cand) > logy) {
        for (int i = 0; i < N; ++i) W(i, j) = cand(i);
        break;
      }
      if (gam < 0.0) gmin = gam; else gmax = gam;
      gam = gmin + (gmax - gmin) * unif_rand();
    }
  }
}

static void update_phi_sep(double& phi, const arma::mat& M, const arma::mat& W,
                           const arma::vec& sig2, const NumericMatrix& locs,
                           const Vecchia& V, Coeff& C, Coeff& Ccand, double nu,
                           double b_phi, double tuning, bool& accepted) {
  accepted = false;
  const double s = std::sqrt(tuning);
  const double post_cur = vecchia_ll_sep(W, M, V, C, sig2);

  double can = -1.0;
  for (int tries = 0; tries < 100000; ++tries) {
    can = phi + s * norm_rand();
    if (can > 0.0 && can < b_phi) break;
    can = -1.0;
  }
  if (can <= 0.0) return;

  vecchia_coeff_fill(locs, V, can, nu, Ccand);
  const double post_can = vecchia_ll_sep(W, M, V, Ccand, sig2);

  const double Zcur = R::pnorm((b_phi - phi) / s, 0, 1, 1, 0) - R::pnorm(-phi / s, 0, 1, 1, 0);
  const double Zcan = R::pnorm((b_phi - can) / s, 0, 1, 1, 0) - R::pnorm(-can / s, 0, 1, 1, 0);
  const double log_r = (post_can - post_cur) + std::log(Zcur) - std::log(Zcan);

  if (std::log(unif_rand()) < log_r) {
    phi = can;
    std::swap(C, Ccand);
    accepted = true;
  }
}

static void update_Sigma_sep(const arma::mat& UW, const arma::mat& M_prior,
                             const arma::mat& cholP, const arma::mat& Mpart,
                             const arma::mat& Rtilde, double a_prior, double b_prior,
                             int N, const std::vector<int>& fam,
                             const std::vector<double>& fixed_sd,
                             arma::vec& sig2) {
  const int q = UW.n_cols;
  for (int j = 0; j < q; ++j) {
    if (fam[j] == 2) {
      sig2(j) = fixed_sd[j] * fixed_sd[j];
      continue;
    }
    const arma::vec rss_p = cholP * M_prior.col(j);
    const arma::vec rss_post = Rtilde * Mpart.col(j);
    double c_val = arma::dot(UW.col(j), UW.col(j)) +
                   arma::dot(rss_p, rss_p) -
                   arma::dot(rss_post, rss_post);
    double shape = a_prior + (double)N / 2.0;
    double rate = b_prior + 0.5 * std::max(c_val, 1e-8);
    sig2(j) = 1.0 / R::rgamma(shape, 1.0 / rate);
  }
}

static arma::mat update_beta_sep(const arma::mat& Mpart, const arma::vec& sig2,
                                 const arma::mat& Rtilde) {
  const int p = Mpart.n_rows, q = Mpart.n_cols;
  arma::mat beta(p, q);
  for (int j = 0; j < q; ++j) {
    arma::vec z(p);
    for (int i = 0; i < p; ++i) z(i) = norm_rand();
    arma::vec tmp = Rtilde * Mpart.col(j) + z * std::sqrt(sig2(j));
    beta.col(j) = Rtilde.t() * tmp;
  }
  return beta;
}

// =============================================================================
// Gibbs sampler for K-fold mixture proposal (separate-covariance model M2)
// =============================================================================
// [[Rcpp::export]]
List gibbs_kfold_sep_cpp(
    const NumericMatrix& Y, const NumericMatrix& X, const NumericMatrix& locs,
    const IntegerMatrix& NNarray, const IntegerVector& fold, int K,
    const List& test_nn_list,
    const CharacterVector& family, const NumericVector& family_par,
    double nu, double phi_init, const NumericMatrix& beta_init,
    const NumericVector& sig2_init, const NumericMatrix& W_init,
    const NumericMatrix& M_prior, const NumericMatrix& V_prior,
    double a_prior, double b_prior,
    double b_phi, double tuning_phi, int niters) {

  const int N = Y.nrow(), q = Y.ncol(), p = X.ncol();

  std::vector<int> fam(q);
  std::vector<double> par(q, 1.0);
  for (int j = 0; j < q; ++j) {
    std::string f = as<std::string>(family[j]);
    if      (f == "Gaussian")          fam[j] = 0;
    else if (f == "Poisson")           fam[j] = 1;
    else if (f == "Binomial")          fam[j] = 2;
    else if (f == "Gamma")             fam[j] = 3;
    else if (f == "Negative-Binomial") fam[j] = 4;
    else stop("Unsupported family: " + f);
    if (family_par.size() == q && !NumericVector::is_na(family_par[j])) par[j] = family_par[j];
  }

  std::vector<int> fold0(N);
  for (int i = 0; i < N; ++i) fold0[i] = fold[i] - 1;

  const arma::mat Ya = as<arma::mat>(Y), Xa = as<arma::mat>(X);
  arma::mat W = as<arma::mat>(W_init);
  arma::mat beta = as<arma::mat>(beta_init);
  arma::vec sig2 = as<arma::vec>(sig2_init);
  const arma::mat Mprior = as<arma::mat>(M_prior);
  const arma::mat precP = arma::inv_sympd(as<arma::mat>(V_prior));
  const arma::mat cholP = arma::chol(precP, "upper");

  std::vector<double> fixed_sd(q, 0.0);
  for (int j = 0; j < q; ++j) fixed_sd[j] = std::sqrt(sig2(j));

  Vecchia V; build_structure(NNarray, V);
  Coeff C, Ccand;
  double phi = phi_init;
  vecchia_coeff_fill(locs, V, phi, nu, C);

  const int ntheta = 1 + p * q + q;
  NumericMatrix samples(niters, ntheta);
  NumericMatrix fold_ll_out(niters, K);
  IntegerVector folds(niters);

  std::vector<double> fold_ll(K);
  compute_all_fold_predictive_loglik_sep(Ya, Xa, W, beta, sig2, locs, test_nn_list,
                                         phi, nu, fam, par, fold_ll);

  std::vector<char> train(N, 1);
  arma::mat M, UW, UX;
  double acc = 0.0;
  const int step_10pct = std::max(1, niters / 10);

  for (int iter = 0; iter < niters; ++iter) {
    // 1. Draw fold index k | theta, W via Gumbel-max trick
    int k = 0;
    double best = R_NegInf;
    for (int l = 0; l < K; ++l) {
      const double v = -fold_ll[l] - std::log(-std::log(unif_rand()));
      if (v > best) { best = v; k = l; }
    }
    folds[iter] = k + 1;
    for (int i = 0; i < N; ++i) train[i] = (fold0[i] != k);

    // 2. Draw theta, W | k ~ p(theta, W | y_{-D_k})
    M = Xa * beta;
    bool accepted = false;
    update_phi_sep(phi, M, W, sig2, locs, V, C, Ccand, nu, b_phi, tuning_phi, accepted);
    acc += accepted;

    vecchia_products(W, Xa, V, C, UW, UX);
    const arma::mat Rtilde = arma::chol(arma::inv_sympd(UX.t() * UX + precP), "upper");
    const arma::mat Mpart = UX.t() * UW + precP * Mprior;

    update_Sigma_sep(UW, Mprior, cholP, Mpart, Rtilde, a_prior, b_prior, N, fam, fixed_sd, sig2);
    beta = update_beta_sep(Mpart, sig2, Rtilde);

    M = Xa * beta;
    update_W_sep(W, M, sig2, V, C, Ya, fam, par, train);

    // 3. Compute predictive log-likelihood in each fold
    compute_all_fold_predictive_loglik_sep(Ya, Xa, W, beta, sig2, locs, test_nn_list,
                                           phi, nu, fam, par, fold_ll);

    for (int l = 0; l < K; ++l) fold_ll_out(iter, l) = fold_ll[l];

    int c = 0;
    samples(iter, c++) = phi;
    for (int j = 0; j < q; ++j) {
      for (int i = 0; i < p; ++i) samples(iter, c++) = beta(i, j);
    }
    for (int j = 0; j < q; ++j) samples(iter, c++) = sig2(j);

    if ((iter & 15) == 0) Rcpp::checkUserInterrupt();
    if ((iter + 1) % step_10pct == 0 || (iter + 1) == niters) {
      Rcout << "Progress: " << (int)((iter + 1) * 100.0 / niters) << "% ("
            << (iter + 1) << " / " << niters << ")\n";
    }
  }

  return List::create(
    _["samples"]     = samples,
    _["theta"]       = samples,
    _["fold_loglik"] = fold_ll_out,
    _["acc"]         = acc / niters,
    _["acc_phi"]     = acc / niters,
    _["folds"]       = folds,
    _["state"]       = List::create(
      _["W"]    = W,
      _["beta"] = beta,
      _["sig2"] = sig2,
      _["phi"]  = phi
    )
  );
}
