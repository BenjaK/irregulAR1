#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

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
arma::vec ar1_sim_cpp(const arma::uword n,
                      const double rho,
                      const double sigma) {
  arma::vec x(n);
  x[0] = R::rnorm(0, 1) * sigma / std::sqrt(1 - std::pow(rho, 2));
  for (auto i = 1; i < n; ++i) {
    x[i] = rho * x[i - 1] + sigma * R::rnorm(0, 1);
  }
  return x;
}

//' Simulate from a stationary Gaussian AR(1) process.
//'
//' Simulate from a stationary Gaussian AR(1) process at \code{n} consecutive
//' time points.
//' @param n The number of timepoints to simulate for.
//' @param mu A vector of length \code{n} with the expected values at each
//'   time point.
//' @param rho A real number strictly less than 1 in absolute value.
//' @param sigma A positive real number.
//' @return A vector of length \code{n} with the process values.
//' @keywords internal
// [[Rcpp::export]]
arma::vec ar1_sim_conditional_cpp(const arma::uvec obs_times,
                                  const arma::uvec pred_times,
                                  const arma::vec obs_x,
                                  const double rho,
                                  const double sigma) {
  auto N = pred_times.n_elem + obs_times.n_elem;
  // Use sort_index on (pred_times, obs_times)
  // Create Q for the sorted vector
  // Draw x_all with this Q; separate into x_pred_star and x_obs_star
  // Create Sigma_po
  // return x_pred_star + Sigma_po * Q * (x_obs - x_obs_star)
  arma::vec x_all(N);
}