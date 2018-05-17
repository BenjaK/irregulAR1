#ifndef SIMULATION_H
#define SIMULATION_H

#include <RcppEigen.h>
#include <cmath>
#include "simulation.h"
#include "matrix.h"
// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd ar1_sim_cpp(const int n,
                      const double rho,
                      const double sigma) {
  Eigen::VectorXd x(n);
  x[0] = R::rnorm(0, 1) * sigma / std::sqrt(1 - std::pow(rho, 2));
  for (auto i = 1; i < n; ++i) {
    x[i] = rho * x[i - 1] + sigma * R::rnorm(0, 1);
  }
  return x;
}

Eigen::VectorXd ar1_sim_irregular_cpp(const Eigen::VectorXi times,
                                      const double rho,
                                      const double sigma) {
  Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(times.size(),
                                                            0.0, 1.0));
  Eigen::SparseMatrix<double> U = chol_tridiag_upper(ar1_prec_irregular(times,
                                                                        rho,
                                                                        sigma));
  return band1_backsolve(U, z);
}

// //' Simulate from a stationary Gaussian AR(1) process.
// //'
// //' Simulate from a stationary Gaussian AR(1) process at \code{n} consecutive
// //' time points.
// //' @param n The number of timepoints to simulate for.
// //' @param mu A vector of length \code{n} with the expected values at each
// //'   time point.
// //' @param rho A real number strictly less than 1 in absolute value.
// //' @param sigma A positive real number.
// //' @return A vector of length \code{n} with the process values.
// //' @keywords internal
// // [[Rcpp::export]]
// Eigen::VectorXd ar1_sim_conditional_cpp(const Eigen::VectorXi pred_times,
//                                   const Eigen::VectorXi obs_times,
//                                   const Eigen::VectorXd obs_x,
//                                   const double rho,
//                                   const double sigma) {
//   auto N = pred_times.size() + obs_times.size();
//   Eigen::VectorXi all_times = arma::join_cols(pred_times, obs_times);
//   Eigen::VectorXi sorted_times = arma::sort_index(all_times);
//   // Use sort_index on (pred_times, obs_times)
//   // Create Q for the sorted vector
//   // Draw x_all with this Q; separate into x_pred_star and x_obs_star
//   // Create Sigma_po
//   // return x_pred_star + Sigma_po * Q * (x_obs - x_obs_star)
//   Eigen::VectorXd x_all(N);
// }

#endif