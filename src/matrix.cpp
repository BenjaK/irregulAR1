#include <RcppEigen.h>
#include <cmath>
#include "matrix.h"
// [[Rcpp::depends(RcppEigen)]]


// Covariance matrices ---------------------------------------------------------

Eigen::MatrixXd ar1_cov_consecutive(const int n,
                                    const double rho,
                                    const double sigma) {
  Eigen::MatrixXd A(n, n);
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

Eigen::MatrixXd ar1_cov_irregular(const Eigen::VectorXi& times,
                                  const double rho,
                                  const double sigma) {
  Eigen::MatrixXd A(times.size(), times.size());
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

Eigen::MatrixXd ar1_cross_cov(const Eigen::VectorXi& times1,
                              const Eigen::VectorXi& times2,
                              const double rho,
                              const double sigma) {
  int n = times1.size();
  int m = times2.size();
  Eigen::MatrixXd A(m, n);

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

Eigen::MatrixXd ar1_cov_chol_irregular(const Eigen::VectorXi& times,
                                       const double rho,
                                       const double sigma) {
  return ar1_cov_irregular(times, rho, sigma).llt().matrixU();
}

// Precision matrices ----------------------------------------------------------

Eigen::SparseMatrix<double> ar1_prec_consecutive(const int n,
                                                 const double rho,
                                                 const double sigma) {
  Eigen::SparseMatrix<double> Q(n, n);
  Q.insert(0, 0) = 1.0;
  Q.insert(n-1, n-1) = 1.0;
  for (int i = 1; i < n; ++i) {
    Q.insert(i, i - 1) = -rho / std::pow(sigma, 2.0);
    Q.insert(i - 1, i) = -rho / std::pow(sigma, 2.0);
    Q.insert(i, i)     = (1.0 + std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  }
  Q.makeCompressed();
  return Q;
}

Eigen::SparseMatrix<double> ar1_prec_irregular(const Eigen::VectorXi& times,
                                               const double rho,
                                               const double sigma) {
  int n = times.size();
  Eigen::SparseMatrix<double> Q(n, n);
  double a = (1.0 - std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  Q.insert(0, 0) = a / (1.0 - std::pow(rho, 2.0 * (times(1) - times(0))));
  Q.insert(n-1, n-1) = a / (1.0 - std::pow(rho, 2.0 * (times(n-1) - times(n-2))));
  double num, den;
  // Fill diagonal
  for (int i = 1; i < n - 1; ++i) {
    num = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i-1)));
    den = (1.0 - std::pow(rho, 2.0 * (times(i) - times(i-1))))
      * (1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i))));
    Q.insert(i, i) = a * num / den;
  }
  // Fill first diagonals above and below main diagonal
  for (int i = 0; i < n - 1; ++i) {
    num = std::pow(rho, times(i+1) - times(i));
    den = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i)));
    Q.insert(i + 1, i) = -a * num / den;
    Q.insert(i, i + 1) = Q.coeff(i + 1, i);
  }
  Q.makeCompressed();
  return Q;
}

Eigen::SparseMatrix<double>
chol_tridiag_upper(const Eigen::SparseMatrix<double>& Q) {
  int n = Q.rows();
  Eigen::SparseMatrix<double> U(n, n);
  int p = 1;
  Eigen::VectorXd v(n);
  for (int j = 0; j < n; ++j) {
    int lambda = std::min(j + p, n - 1);
    // Eigen::VectorXd Q_sub(lambda - j + 1);
    for (auto k = j, i = 0; k <= lambda; ++k, ++i) {
      v(k) = Q.coeff(k, j);
      // Q_sub(i) = Q.coeff(k, j);
    }
    // v.subvec(j, lambda) = Q_sub;
    for (int k = std::max(0, j - p); k < j; ++k) {
      int i = std::min(k + p, n - 1);
      for (int h = 0; h <= i; ++h) {
        v(h) -= U.coeff(k, h) * U.coeff(k, j);
      }
      // v.subvec(j, i) = v.subvec(j, i) - U(k, arma::span(j, i)).t() * U(k, j);
    }
    for (int k = j; k <= lambda; ++k) {
      U.coeffRef(j, k) = v(k) / std::sqrt(v(j));
    }
    // U(j, arma::span(j, lambda)) = v.subvec(j, lambda).t() / std::sqrt(v(j));
  }
  Q.makeCompressed();
  return U;
}

Eigen::VectorXd band1_backsolve(const Eigen::SparseMatrix<double>& U,
                                const Eigen::VectorXd& z) {
  int m = z.size() - 1;
  Eigen::VectorXd v(m);
  v(m) = z(m) / U.coeff(m, m);
  for (int i = m - 2; i >= 0; --i) {
    v(i) = (z(i) - U.coeff(i, i + 1)) / U.coeff(i, i);
  }
  return v;
}

