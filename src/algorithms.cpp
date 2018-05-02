// #include <RcppEigen.h>
// #include <cmath>
// using namespace Rcpp;
// // [[Rcpp::depends(RcppEigen)]]
//
// // Covariance matrices ---------------------------------------------------------
//
// //' Covariance matrix of vector x observed at n consecutive times.
// //' @export
// // [[Rcpp::export]]
// Eigen::MatrixXd ar1_cov_consecutive(const unsigned int n,
//                                     const double rho,
//                                     const double sigma) {
//   Eigen::MatrixXd A(n, n);
//   double t;
//   for (int j = 0; j < n; ++j) {
//     for (int i = j + 1; i < n; ++i) {
//       t = static_cast<double>(i - j);
//       A(i, j) = std::pow(rho, t) * pow(sigma, 2) / (1 - std::pow(rho, 2));
//       A(j, i) = A(i, j);
//     }
//     A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
//   }
//   return A;
// }
//
// //' Covariance matrix of vector x observed at (irregularly spaced) times times.
// //' @export
// // [[Rcpp::export]]
// Eigen::MatrixXd ar1_cov_irregular(const Eigen::VectorXi& times,
//                                   const double rho,
//                                   const double sigma) {
//   Eigen::MatrixXd A(times.size(), times.size());
//   double t1, t2;
//   for (int j = 0; j < times.size(); ++j) {
//     t1 = static_cast<double>(times(j));
//     for (int i = j + 1; i < times.size(); ++i) {
//       t2 = static_cast<double>(times(i));
//       A(i, j) = std::pow(rho, std::abs(t1 - t2))
//               * pow(sigma, 2) / (1 - std::pow(rho, 2));
//       A(j, i) = A(i, j);
//     }
//     A(j, j) = pow(sigma, 2) / (1 - std::pow(rho, 2));
//   }
//   return A;
// }