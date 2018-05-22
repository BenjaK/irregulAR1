#include <RcppArmadillo.h>
#include <cmath>
#include "matrix.h"
// [[Rcpp::depends(RcppArmadillo)]]


// Covariance matrices ---------------------------------------------------------

arma::mat ar1_cov_consecutive(const int n,
                              const double rho,
                              const double sigma) {
  arma::mat A(n, n);
  double t;
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      t = static_cast<double>(i - j);
      A(i, j) = std::pow(rho, t) * pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

arma::mat ar1_cov_irregular(const arma::uvec& times,
                            const double rho,
                            const double sigma) {
  arma::mat A(times.size(), times.size());
  double t1, t2;
  for (int j = 0; j < times.size(); ++j) {
    t1 = static_cast<double>(times(j));
    for (int i = j + 1; i < times.size(); ++i) {
      t2 = static_cast<double>(times(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
                pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

arma::mat ar1_cross_cov(const arma::uvec& times1,
                        const arma::uvec& times2,
                        const double rho,
                        const double sigma) {
  int n = times1.size();
  int m = times2.size();
  arma::mat A(m, n);

  double t1, t2;
  for (int j = 0; j < n; ++j) {
    t1 = static_cast<double>(times1(j));
    for (int i = 0; i < m; ++i) {
      t2 = static_cast<double>(times2(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
        pow(sigma, 2) / (1 - std::pow(rho, 2));
    }
  }
  return A;
}

arma::mat ar1_cov_chol_irregular(const arma::uvec& times,
                                 const double rho,
                                 const double sigma) {
  return arma::chol(ar1_cov_irregular(times, rho, sigma));
}

// Precision matrices ----------------------------------------------------------

arma::sp_mat ar1_prec_consecutive(const int n,
                                  const double rho,
                                  const double sigma) {
  arma::sp_mat Q(n, n);
  Q(0, 0) = 1.0 / std::pow(sigma, 2.0);
  for (int i = 1; i < n; ++i) {
    Q(i, i - 1) = -rho / std::pow(sigma, 2.0);
    Q(i - 1, i) = -rho / std::pow(sigma, 2.0);
    Q(i, i)     = (1.0 + std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  }
  Q(n-1, n-1) = 1.0 / std::pow(sigma, 2.0);
  return Q;
}

arma::sp_mat ar1_prec_irregular(const arma::uvec& times,
                                const double rho,
                                const double sigma) {
  int n = times.size();
  arma::sp_mat Q(n, n);
  double a = (1.0 - std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  Q(0, 0) = a / (1.0 - std::pow(rho, 2.0 * (times(1) - times(0))));
  Q(n-1, n-1) = a / (1.0 - std::pow(rho, 2.0 * (times(n-1) - times(n-2))));
  double num, den;
  // Fill diagonal
  for (int i = 1; i < n - 1; ++i) {
    num = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i-1)));
    den = (1.0 - std::pow(rho, 2.0 * (times(i) - times(i-1))))
      * (1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i))));
    Q(i, i) = a * num / den;
  }
  // Fill first diagonals above and below main diagonal
  for (int i = 0; i < n - 1; ++i) {
    num = std::pow(rho, times(i+1) - times(i));
    den = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i)));
    Q(i + 1, i) = -a * num / den;
    Q(i, i + 1) = Q(i + 1, i);
  }
  return Q;
}

arma::sp_mat chol_tridiag_upper(const arma::sp_mat& Q) {
  int n = Q.n_rows;
  arma::sp_mat U(n, n);
  int p = 1;
  arma::vec v(n);
  for (int j = 0; j < n; ++j) {
    int lambda = std::min(j + p, n - 1);
    for (auto k = j, i = 0; k <= lambda; ++k, ++i) {
      v(k) = Q(k, j);
    }
    for (int k = std::max(0, j - p); k < j; ++k) {
      int i = std::min(k + p, n - 1);
      for (int h = 0; h <= i; ++h) {
        v(h) -= U(k, h) * U(k, j);
      }
    }
    U(j, arma::span(j, lambda)) = v.subvec(j, lambda).t() / std::sqrt(v(j));
  }
  return U;
}

arma::vec band1_backsolve_cpp(const arma::sp_mat& U,
                              const arma::vec& z) {
  int m = z.size();
  arma::vec v(m);
  v(m-1) = z(m-1) / U(m-1, m-1);
  for (int i = m - 2; i >= 0; --i) {
    v(i) = (z(i) - U(i, i + 1) * v(i+1)) / U(i, i);
  }
  return v;
}

