#include <Rcpp.h>
using namespace Rcpp;

// Numerically stable log-sum-exp
double log_sum_exp(const NumericVector& x) {
  double m = Rcpp::max(x);
  
  double s = 0.0;
  for (int i = 0; i < x.size(); ++i)
    s += std::exp(x[i] - m);
  
  return m + std::log(s);
}


// Compute potential U and gradient simultaneously
// for the LOOCV mixture proposal.
//
// [[Rcpp::export]]
List eval_loocv(
    const NumericMatrix& X,
    const NumericVector& y,
    const NumericVector& theta,
    double sigma2,
    double tau2
) {
  
  int n = X.nrow();
  int d = X.ncol();
  
  NumericVector resid(n);
  NumericVector ll(n);
  
  // --------------------------------------------------
  // First pass:
  // compute residuals and observation log likelihoods
  // --------------------------------------------------
  
  for (int i = 0; i < n; ++i) {
    
    double mu = 0.0;
    
    for (int j = 0; j < d; ++j)
      mu += X(i,j) * theta[j];
    
    resid[i] = y[i] - mu;
    
    ll[i] =
      -0.5 * std::log(2.0 * M_PI * sigma2)
      -0.5 * resid[i] * resid[i] / sigma2;
  }
  
  // --------------------------------------------------
  // log mixture term
  // --------------------------------------------------
  
  NumericVector neg_ll(n);
  
  for (int i = 0; i < n; ++i)
    neg_ll[i] = -ll[i];
  
  double lse = log_sum_exp(neg_ll);
  
  // --------------------------------------------------
  // Gradient
  // --------------------------------------------------
  
  NumericVector grad(d);
  
  // prior contribution
  double log_prior = 0.0;
  
  for (int j = 0; j < d; ++j) {
    log_prior +=
      -0.5 * std::log(2.0 * M_PI * tau2)
      -0.5 * theta[j] * theta[j] / tau2;
      
      grad[j] = -theta[j] / tau2;
  }
  
  // ordinary likelihood contribution
  for (int i = 0; i < n; ++i) {
    
    double w = std::exp(-ll[i] - lse);
    double r = resid[i] / sigma2;
    
    for (int j = 0; j < d; ++j) {
      grad[j] += X(i,j) * r * (1.0 - w);
    }
  }
  
  // log q
  double log_q =
    log_prior;
  
  for (int i = 0; i < n; ++i)
    log_q += ll[i];
  
  log_q += lse;
  
  // U = -log q
  double U = -log_q;
  
  // grad_U = -grad_log_q
  grad = -grad;
  
  return List::create(
    _["U"] = U,
    _["grad"] = grad
  );
}


// [[Rcpp::export]]
List hmc_loocv_cpp(
    const NumericMatrix& X,
    const NumericVector& y,
    const NumericVector& theta_init,
    double sigma2,
    double tau2,
    double epsilon,
    int L,
    int niters
) {
  
  int d = X.ncol();
  
  NumericMatrix samples(niters, d);
  
  NumericVector theta = clone(theta_init);
  
  for (int j = 0; j < d; ++j)
    samples(0,j) = theta[j];
  
  int acc = 0;
  
  // Evaluate current state ONCE
  List current = eval_loocv(
    X, y, theta, sigma2, tau2
  );
  
  double U_current = current["U"];
  NumericVector grad_current = current["grad"];
  
  for (int iter = 1; iter < niters; ++iter) {
    
    // Fresh momentum
    NumericVector p(d);
    
    for (int j = 0; j < d; ++j)
      p[j] = R::rnorm(0.0, 1.0);
    
    NumericVector theta_prop = clone(theta);
    NumericVector p_prop = clone(p);
    NumericVector grad = clone(grad_current);
    
    // Initial half-step
    for (int j = 0; j < d; ++j)
      p_prop[j] -= 0.5 * epsilon * grad[j];
    
    double U_prop = NA_REAL;
    
    // Leapfrog
    for (int l = 0; l < L; ++l) {
      
      for (int j = 0; j < d; ++j)
        theta_prop[j] += epsilon * p_prop[j];
      
      List val = eval_loocv(
        X, y, theta_prop, sigma2, tau2
      );
      
      U_prop = val["U"];
      grad = val["grad"];
      
      if (l != L - 1) {
        for (int j = 0; j < d; ++j)
          p_prop[j] -= epsilon * grad[j];
      }
    }
    
    // Final half-step
    for (int j = 0; j < d; ++j)
      p_prop[j] -= 0.5 * epsilon * grad[j];
    
    // Hamiltonians
    double K_current = 0.0;
    double K_prop = 0.0;
    
    for (int j = 0; j < d; ++j) {
      K_current += 0.5 * p[j] * p[j];
      K_prop += 0.5 * p_prop[j] * p_prop[j];
    }
    
    double log_alpha =
      -(U_prop + K_prop)
      +(U_current + K_current);
      
      if (std::log(R::runif(0.0, 1.0)) < log_alpha) {
        
        theta = theta_prop;
        grad_current = grad;
        U_current = U_prop;
        
        acc++;
        
      } else {
        // theta, grad_current, U_current unchanged
      }
      
      for (int j = 0; j < d; ++j)
        samples(iter,j) = theta[j];
  }
  
  return List::create(
    _["samples"] = samples,
    _["acc"] = (double) acc / (niters - 1)
  );
}
