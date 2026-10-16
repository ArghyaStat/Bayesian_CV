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

static double vecchia_ll(const arma::mat& W, const arma::mat& M, const Vecchia& V,
                         const Coeff& C, const arma::mat& cholV) {
  const int N = W.n_rows, q = W.n_cols;
  const arma::mat E = W - M;
  double logdetV = 0.0;
  for (int j = 0; j < q; ++j) logdetV += 2.0 * std::log(cholV(j, j));

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
    for (int j = 0; j < q; ++j) {
      double t = res[j];
      for (int l = 0; l < j; ++l) t -= cholV(l, j) * z[l];
      z[j] = t / cholV(j, j);
    }
    double quad = 0.0;
    for (int j = 0; j < q; ++j) quad += z[j] * z[j];
    const double ri = C.r[i];
    ll += -0.5 * (q * std::log(2.0 * M_PI * ri) + logdetV + quad / ri);
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

// Out-of-sample Vecchia predictive log-likelihood for all folds
static void compute_all_fold_predictive_loglik(
    const arma::mat& Y, const arma::mat& X, const arma::mat& W,
    const arma::mat& beta, const arma::mat& cholSigma,
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

      arma::vec z(q);
      for (int j = 0; j < q; ++j) z(j) = norm_rand();
      arma::rowvec w_draw = mean_i + std::sqrt(r_i) * (z.t() * cholSigma);

      for (int j = 0; j < q; ++j) {
        ll_fold += loglik1(fam[j], par[j], Y(site_i, j), w_draw(j));
      }
    }
    out[k] = ll_fold;
  }
}

static void update_W(arma::mat& W, const arma::mat& M, const arma::mat& Sigma,
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
    double Scond = Sigma(j, j);
    arma::vec cvec;
    arma::uvec minus;
    if (q > 1) {
      minus.set_size(q - 1);
      int c = 0;
      for (int l = 0; l < q; ++l) if (l != j) minus(c++) = l;
      arma::mat Smm = Sigma.submat(minus, minus);
      arma::vec Smj = Sigma.submat(minus, arma::uvec({(arma::uword)j}));
      cvec = arma::solve(Smm, Smj);
      Scond = Sigma(j, j) - arma::dot(Smj, cvec);
    }
    for (int i = 0; i < N; ++i) {
      double m = M(i, j);
      if (q > 1) {
        for (int c = 0; c < q - 1; ++c) {
          m += (W(i, minus(c)) - M(i, minus(c))) * cvec(c);
        }
      }
      mu(i) = m;
    }
    const double sc = std::sqrt(std::max(Scond, 0.0));

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

static void update_phi(double& phi, const arma::mat& M, const arma::mat& W,
                       const arma::mat& cholSigma, const NumericMatrix& locs,
                       const Vecchia& V, Coeff& C, Coeff& Ccand, double nu,
                       double b_phi, double tuning, bool& accepted) {
  accepted = false;
  const double s = std::sqrt(tuning);
  const double post_cur = vecchia_ll(W, M, V, C, cholSigma);

  double can = -1.0;
  for (int tries = 0; tries < 100000; ++tries) {
    can = phi + s * norm_rand();
    if (can > 0.0 && can < b_phi) break;
    can = -1.0;
  }
  if (can <= 0.0) return;

  vecchia_coeff_fill(locs, V, can, nu, Ccand);
  const double post_can = vecchia_ll(W, M, V, Ccand, cholSigma);

  const double Zcur = R::pnorm((b_phi - phi) / s, 0, 1, 1, 0) - R::pnorm(-phi / s, 0, 1, 1, 0);
  const double Zcan = R::pnorm((b_phi - can) / s, 0, 1, 1, 0) - R::pnorm(-can / s, 0, 1, 1, 0);
  const double log_r = (post_can - post_cur) + std::log(Zcur) - std::log(Zcan);

  if (std::log(unif_rand()) < log_r) {
    phi = can;
    std::swap(C, Ccand);
    accepted = true;
  }
}

static arma::mat riwish_cpp(double df, const arma::mat& S) {
  const int q = S.n_rows;
  arma::mat Sinv = arma::inv_sympd(S);
  arma::mat R = arma::chol(Sinv, "upper");
  arma::mat A(q, q, arma::fill::zeros);
  for (int i = 0; i < q; ++i) {
    A(i, i) = std::sqrt(R::rchisq(df - i));
    for (int j = i + 1; j < q; ++j) {
      A(i, j) = norm_rand();
    }
  }
  arma::mat AR = A * R;
  arma::mat W = AR.t() * AR;
  return arma::inv_sympd(W);
}

static void update_Sigma(const arma::mat& UW, const arma::mat& M_prior,
                         const arma::mat& cholP, const arma::mat& Mpart,
                         const arma::mat& Rtilde, const arma::mat& S_prior,
                         double df_prior, int N, const std::vector<int>& fam,
                         const std::vector<double>& fixed_sd,
                         arma::mat& Sigma, arma::mat& cholSigma) {
  const int q = UW.n_cols;
  const arma::mat rssW = UW.t() * UW;
  const arma::mat a = cholP * M_prior;
  const arma::mat b = Rtilde * Mpart;
  arma::mat S_tilde = S_prior + rssW + a.t() * a - b.t() * b;

  const arma::mat Sstar = riwish_cpp(df_prior + N, S_tilde);
  arma::vec dstar = arma::sqrt(Sstar.diag());
  arma::mat Rstar = Sstar.each_col() / dstar;
  Rstar.each_row() /= dstar.t();

  arma::vec dnew = dstar;
  for (int j = 0; j < q; ++j) {
    if (fam[j] == 2) dnew(j) = fixed_sd[j];
  }

  Sigma = Rstar;
  Sigma.each_col() %= dnew;
  Sigma.each_row() %= dnew.t();
  cholSigma = arma::chol(Sigma, "upper");
}

static arma::mat update_beta(const arma::mat& Mpart, const arma::mat& cholSigma,
                             const arma::mat& Rtilde) {
  const int p = Mpart.n_rows, q = Mpart.n_cols;
  arma::mat Z(p, q);
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < p; ++i) Z(i, j) = norm_rand();
  }
  arma::mat tmp = Rtilde * Mpart + Z * cholSigma;
  return Rtilde.t() * tmp;
}

// [[Rcpp::export]]
List gibbs_kfold_joint_cpp(
    const NumericMatrix& Y, const NumericMatrix& X, const NumericMatrix& locs,
    const IntegerMatrix& NNarray, const IntegerVector& fold, int K,
    const List& test_nn_list,
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
  for (int i = 0; i < N; ++i) {
    fold0[i] = fold[i] - 1;
    if (fold0[i] < 0 || fold0[i] >= K) stop("fold must take values in 1..K.");
  }

  const arma::mat Ya = as<arma::mat>(Y), Xa = as<arma::mat>(X);
  arma::mat W = as<arma::mat>(W_init);
  arma::mat beta = as<arma::mat>(beta_init);
  arma::mat Sigma = as<arma::mat>(Sigma_init);
  arma::mat cholSigma = arma::chol(Sigma, "upper");
  const arma::mat Mprior = as<arma::mat>(M_prior);
  const arma::mat Sprior = as<arma::mat>(S_prior);

  const arma::mat Vprior = as<arma::mat>(V_prior);
  const arma::mat precP = arma::inv_sympd(Vprior);
  const arma::mat cholP = arma::chol(precP, "upper");

  std::vector<double> fixed_sd(q, 0.0);
  for (int j = 0; j < q; ++j) fixed_sd[j] = std::sqrt(Sigma(j, j));

  Vecchia V; build_structure(NNarray, V);
  Coeff C, Ccand;
  double phi = phi_init;
  vecchia_coeff_fill(locs, V, phi, nu, C);

  const int ntheta = 1 + p * q + q * (q + 1) / 2;
  NumericMatrix fold_ll_out(niters, K);
  NumericMatrix theta_out(niters, ntheta);
  IntegerVector k_out(niters);

  std::vector<double> fold_ll(K);
  compute_all_fold_predictive_loglik(Ya, Xa, W, beta, cholSigma, locs, test_nn_list,
                                     phi, nu, fam, par, fold_ll);

  std::vector<char> train(N, 1);
  arma::mat M, UW, UX;

  int acc_main = 0;
  int k = 0;
  const int step_10pct = std::max(1, niters / 10);
  auto t_start = std::chrono::steady_clock::now();

  for (int iter = 0; iter < niters; ++iter) {

    // 1. Draw fold index kappa | theta via Gumbel-max trick
    double best = R_NegInf;
    for (int l = 0; l < K; ++l) {
      const double v = -fold_ll[l] - std::log(-std::log(unif_rand()));
      if (v > best) {
        best = v;
        k = l;
      }
    }
    for (int i = 0; i < N; ++i) train[i] = (fold0[i] != k);

    // 2. Draw theta | kappa ~ p(theta | y_{-F_kappa})
    M = Xa * beta;

    bool acc = false;
    update_phi(phi, M, W, cholSigma, locs, V, C, Ccand, nu, b_phi, tuning_phi, acc);
    acc_main += acc;

    vecchia_products(W, Xa, V, C, UW, UX);
    const arma::mat prec_tilde = UX.t() * UX + precP;
    const arma::mat V_tilde = arma::inv_sympd(prec_tilde);
    const arma::mat Rtilde = arma::chol(V_tilde, "upper");
    const arma::mat Mpart = UX.t() * UW + precP * Mprior;

    update_Sigma(UW, Mprior, cholP, Mpart, Rtilde, Sprior, df_prior, N, fam, fixed_sd,
                 Sigma, cholSigma);
    beta = update_beta(Mpart, cholSigma, Rtilde);

    M = Xa * beta;
    update_W(W, M, Sigma, V, C, Ya, fam, par, train);

    // 3. Compute predictive log-likelihood in each fold at the new theta
    compute_all_fold_predictive_loglik(Ya, Xa, W, beta, cholSigma, locs, test_nn_list,
                                       phi, nu, fam, par, fold_ll);

    for (int l = 0; l < K; ++l) fold_ll_out(iter, l) = fold_ll[l];
    k_out[iter] = k + 1;
    int c = 0;
    theta_out(iter, c++) = phi;
    for (int j = 0; j < q; ++j) {
      for (int i = 0; i < p; ++i) theta_out(iter, c++) = beta(i, j);
    }
    for (int j = 0; j < q; ++j) {
      for (int i = 0; i <= j; ++i) theta_out(iter, c++) = Sigma(i, j);
    }

    if ((iter & 15) == 0) Rcpp::checkUserInterrupt();

    if ((iter + 1) % step_10pct == 0 || (iter + 1) == niters) {
      auto t_now = std::chrono::steady_clock::now();
      double elapsed_mins = std::chrono::duration<double>(t_now - t_start).count() / 60.0;
      int pct = (int)((iter + 1) * 100.0 / niters);
      Rcout << "Progress: " << pct << "% (" << (iter + 1) << " / " << niters
            << ") - Elapsed: " << std::setprecision(2) << std::fixed << elapsed_mins << " mins\n";
    }
  }

  double rate_main = niters > 0 ? (double)acc_main / niters : NA_REAL;

  return List::create(
    _["samples"]     = theta_out,
    _["theta"]       = theta_out,
    _["fold_loglik"] = fold_ll_out,
    _["acc"]         = rate_main,
    _["acc.phi"]     = rate_main,
    _["acc_phi"]     = rate_main,
    _["folds"]       = k_out,
    _["state"]       = List::create(
      _["W"]     = W,
      _["beta"]  = beta,
      _["Sigma"] = Sigma,
      _["phi"]   = phi
    )
  );
}

