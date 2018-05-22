#' Evaluate the log-density of a stationary Gaussian AR(1) process.
#'
#' Evaluate the log-density of a stationary Gaussian AR(1) process, observed at
#' times \code{times} taking values \code{x}.
#' @param x A vector of observed values.
#' @param mu A vector of expected values.
#' @param times A vector of the time points of observation.
#' @param rho A real number strictly less than 1 in absolute value.
#' @param sigma A positive real number.
#' @return A scalar, the log density.
#' @export
ar1_lpdf <- function(x, times, rho, sigma, mu = 0) {
  if (length(mu) != 1 || length(mu) != length(x) || length(x) != length(times))
    stop("Lengths of x, times, and mu (if not a single value) must match.")
  if (abs(rho) >= 1) stop("rho must be less than 1 in magnitude.")
  if (sigma < 0) stop("sigma must be positive")

  return(ar1_lpdf_cpp(x, mu, times, rho, sigma))
}