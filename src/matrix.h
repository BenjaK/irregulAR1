#ifndef MATRIX_H
#define MATRIX_H

#include <RcppEigen.h>
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]

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
Eigen::MatrixXd ar1_cov_consecutive(const int n,
                                    const double rho,
                                    const double sigma);

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
Eigen::MatrixXd ar1_cov_irregular(const Eigen::VectorXi& times,
                                  const double rho,
                                  const double sigma);

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
Eigen::MatrixXd ar1_cross_cov(const Eigen::VectorXi& times1,
                              const Eigen::VectorXi& times2,
                              const double rho,
                              const double sigma);

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
Eigen::MatrixXd ar1_cov_chol_irregular(const Eigen::VectorXi& times,
                                       const double rho,
                                       const double sigma);

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
Eigen::SparseMatrix<double> ar1_prec_consecutive(const int n,
                                                 const double rho,
                                                 const double sigma);

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
Eigen::SparseMatrix<double> ar1_prec_irregular(const Eigen::VectorXi& times,
                                               const double rho,
                                               const double sigma);

//' Upper Cholesky decomposition of a tridiagonal matrix.
//'
//' Creates the lower Cholesky decomposition of a tridiagonal matrix. The
//' decomposition will be a sparse lower triangular matrix with non-zero
//' elements only on the main diagonal and the diagonal below it.
//' @param Q A square tridiagonal matrix.
//' @return A sparse square matrix with the same size as the input matrix.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double>
chol_tridiag_upper(const Eigen::SparseMatrix<double>& Q);


//' Backsolve with band 1 upper Cholesky.
//'
//' Backsolve with band 1 upper Cholesky.
//' @param Q A square tridiagonal matrix.
//' @return A sparse square matrix with the same size as the input matrix.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd band1_backsolve(const Eigen::SparseMatrix<double>& U,
                                const Eigen::VectorXd& z);

#endif
