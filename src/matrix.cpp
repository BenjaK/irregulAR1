#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// Covariance matrices ---------------------------------------------------------

//' Covariance matrix for a stationary Gaussian AR(1) process, observed at
//' consecutive timepoints.
//'
//' Creates the covariance matrix of an AR(1) process with parameters \code{rho}
//' and \code{sigma}, observed at \code{n} consecutive time points. The process
//' is assumed to be in stationarity and have Gaussian errors.
//' @param n An integer greater than or equal to 1.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{n} rows and \code{n} columns.
//' @export
// [[Rcpp::export]]
arma::mat ar1_cov_consecutive(const arma::uword n,
                              const double rho,
                              const double sigma) {
  arma::mat A(n, n);
  double t;
  for (arma::sword j = 0; j < n; ++j) {
    for (arma::sword i = j + 1; i < n; ++i) {
      t = static_cast<double>(i - j);
      A(i, j) = std::pow(rho, t) * pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

//' Covariance matrix for a stationary Gaussian AR(1) process, observed at
//' irregularly spaced time points.
//'
//' Creates the covariance matrix of an AR(1) process with parameters \code{rho}
//' and \code{sigma}, observed at the time points in the vector \code{times}.
//' The process is assumed to be in stationarity and have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
// [[Rcpp::export]]
arma::mat ar1_cov_irregular(const arma::uvec& times,
                            const double rho,
                            const double sigma) {
  arma::mat A(times.n_elem, times.n_elem);
  double t1, t2;
  for (arma::sword j = 0; j < times.n_elem; ++j) {
    t1 = static_cast<double>(times(j));
    for (arma::sword i = j + 1; i < times.n_elem; ++i) {
      t2 = static_cast<double>(times(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
                pow(sigma, 2) / (1 - std::pow(rho, 2));
      A(j, i) = A(i, j);
    }
    A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
  }
  return A;
}

//' Cross-covariance matrix of a stationary Gaussian AR(1) process.
//'
//' Creates the cross-covariance matrix of an AR(1) process with parameters
//' \code{rho} and \code{sigma}, observed at (positive) integer times
//' \code{times1} and \code{times2}, which may be irregularly spaced. The
//' process is assumed to be in stationarity and have Gaussian errors.
//' @param times1 An vector of positive integers, preferably ordered.
//' @param times2 An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{length(times2)} rows and \code{length(times1)}
//'   columns.
//' @export
// [[Rcpp::export]]
arma::mat ar1_cross_cov(const arma::uvec& times1,
                        const arma::uvec& times2,
                        const double rho,
                        const double sigma) {
  arma::uword n = times1.n_elem;
  arma::uword m = times2.n_elem;
  arma::mat A(m, n);

  double t1, t2;
  for (arma::sword j = 0; j < n; ++j) {
    t1 = static_cast<double>(times1(j));
    for (arma::sword i = 0; i < m; ++i) {
      t2 = static_cast<double>(times2(i));
      A(i, j) = std::pow(rho, std::abs(t1 - t2)) *
        pow(sigma, 2) / (1 - std::pow(rho, 2));
    }
  }
  return A;
}

//' Upper triangular Cholesky decomposition for a stationary Gaussian AR(1)
//' process, observed at irregularly spaced time points.
//'
//' Creates the upper triangular Cholesky decomposition matrix of an AR(1)
//' process with parameters \code{rho} and \code{sigma}, observed at the time
//' points in the vector \code{times}. The process is assumed to be in
//' stationarity and have Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
// [[Rcpp::export]]
arma::mat ar1_cov_chol_irregular(const arma::uvec& times,
                                 const double rho,
                                 const double sigma) {
  return arma::chol(ar1_cov_irregular(times, rho, sigma));
}

// Precision matrices ----------------------------------------------------------

//' Sparse precision matrix for a stationary Gaussian AR(1) process, observed at
//' consecutive timepoints.
//'
//' Creates the precision (inverse covariance) matrix of an AR(1) process with
//' parameters \code{rho} and \code{sigma}, observed at \code{n} consecutive
//' time points. The process is assumed to be in stationarity and have Gaussian
//' errors. The matrix is a tridiagonal band matrix and thus sparse.
//' @param n An integer greater than or equal to 1.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A matrix with \code{n} rows and \code{n} columns.
//' @export
// [[Rcpp::export]]
arma::sp_mat ar1_prec_consecutive(const arma::uword n,
                                  const double rho,
                                  const double sigma) {
  arma::sp_mat Q(n, n);
  Q(0, 0) = 1.0;
  Q(n-1, n-1) = 1.0;
  for (arma::uword i = 1; i < n; ++i) {
    Q(i, i - 1) = -rho / std::pow(sigma, 2.0);
    Q(i - 1, i) = -rho / std::pow(sigma, 2.0);
    Q(i, i)     = (1.0 + std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  }
  return Q;
}

//' Precision matrix for a stationary Gaussian AR(1) process, observed at
//' irregularly spaced time points.
//'
//' Creates the precision (inverse covariance) matrix of an AR(1) process with
//' parameters \code{rho} and \code{sigma}, observed at the time points in the
//' vector \code{times}. The process is assumed to be in stationarity and have
//' Gaussian errors.
//' @param times An vector of positive integers, preferably ordered.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A square matrix with \code{length(times)} rows.
//' @export
// [[Rcpp::export]]
arma::sp_mat ar1_prec_irregular(const arma::uvec& times,
                                const double rho,
                                const double sigma) {
  arma::uword n = times.n_elem;
  arma::sp_mat Q(n, n);
  double a = (1.0 - std::pow(rho, 2.0)) / std::pow(sigma, 2.0);
  Q(0, 0) = a / (1.0 - std::pow(rho, 2.0 * (times(1) - times(0))));
  Q(n-1, n-1) = a / (1.0 - std::pow(rho, 2.0 * (times(n-1) - times(n-2))));
  double num, den;
  // Fill diagonal
  for (arma::uword i = 1; i < n - 1; ++i) {
    num = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i-1)));
    den = (1.0 - std::pow(rho, 2.0 * (times(i) - times(i-1))))
      * (1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i))));
    Q(i, i) = a * num / den;
  }
  // Fill first diagonals above and below main diagonal
  for (arma::uword i = 0; i < n - 1; ++i) {
    num = std::pow(rho, times(i+1) - times(i));
    den = 1.0 - std::pow(rho, 2.0 * (times(i+1) - times(i)));
    Q(i + 1, i) = -a * num / den;
    Q(i, i + 1) = Q(i + 1, i);
  }
  return Q;
}

//' Lower Cholesky decomposition of a tridiagonal matrix.
//'
//' Creates the lower Cholesky decomposition of a tridiagonal matrix. The
//' decomposition will be a sparse lower triangular matrix with non-zero
//' elements only on the main diagonal and the diagonal below it.
//' @param Q A square tridiagonal matrix.
//' @return A sparse square matrix with the same size as the input matrix.
//' @export
// [[Rcpp::export]]
arma::sp_mat chol_tridiag_lower(const arma::mat& Q) {
  int n = Q.n_rows;
  arma::sp_mat L(n, n);
  arma::vec v = arma::zeros<arma::vec>(n);
  for (int j = 0; j < n; ++j) {
    int lambda = std::min(j + 1, n - 1);
    v.subvec(j, lambda) = Q(arma::span(j, lambda), j);
    if (j > 0) {
      v(j) -= std::pow(L(j, j - 1), 2.0);
    }
    L(arma::span(j, lambda), j) = v.subvec(j, lambda) / std::sqrt(v(j));
  }
  return L;
}

//' Upper Cholesky decomposition of a tridiagonal matrix.
//'
//' Creates the lower Cholesky decomposition of a tridiagonal matrix. The
//' decomposition will be a sparse lower triangular matrix with non-zero
//' elements only on the main diagonal and the diagonal below it.
//' @param Q A square tridiagonal matrix.
//' @return A sparse square matrix with the same size as the input matrix.
//' @export
// [[Rcpp::export]]
arma::sp_mat chol_tridiag_upper(const arma::mat& Q) {
  int n = Q.n_rows;
  arma::sp_mat U(n, n);
  int p = 1;
  arma::vec v(n);
  for (int j = 0; j < n; ++j) {
    int lambda = std::min(j + p, n - 1);
    v.subvec(j, lambda) = Q(arma::span(j, lambda), j);
    for (int k = std::max(0, j - p); k < j; ++k) {
      int i = std::min(k + p, n - 1);
      v.subvec(j, i) = v.subvec(j, i) - U(k, arma::span(j, i)).t() * U(k, j);
    }
    U(j, arma::span(j, lambda)) = v.subvec(j, lambda).t() / std::sqrt(v(j));
  }
  return U;
}

