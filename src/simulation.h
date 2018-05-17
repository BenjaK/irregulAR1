#ifndef SIMULATION_H
#define SIMULATION_H

#include <RcppEigen.h>
#include <cmath>
#include "matrix.h"
// [[Rcpp::depends(RcppEigen)]]


//' Simulate from a stationary Gaussian AR(1) process.
//'
//' Simulate from a stationary Gaussian AR(1) process at \code{n} consecutive
//' time points.
//' @param n The number of timepoints to simulate for.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A vector of length \code{n} with the process values.
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd ar1_sim_cpp(const int n,
                            const double rho,
                            const double sigma);

//' Simulate from a stationary Gaussian AR(1) process at irregular times.
//'
//' Simulate from a stationary Gaussian AR(1) process at irregular times.
//' @param n The number of timepoints to simulate for.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A vector of length \code{n} with the process values.
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd ar1_sim_irregular_cpp(const Eigen::VectorXi times,
                                      const double rho,
                                      const double sigma);

#endif