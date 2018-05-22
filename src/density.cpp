#include <cmath>
#include "density.h"
#include "matrix.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double ar1_lpdf_cpp(const arma::vec& x,
                    const arma::vec& mu,
                    const arma::uvec& times,
                    const double rho,
                    const double sigma) {
  arma::sp_mat Q = ar1_prec_irregular(times, rho, sigma);
  arma::sp_mat U = chol_tridiag_upper(Q);

  arma::vec z = U * (x - mu);
  double q = arma::dot(z, z);
  return - x.n_elem / 2 * std::log(2 * arma::datum::pi)
         + arma::accu(arma::log(arma::vec(U.diag())))
         - q / 2;
}
