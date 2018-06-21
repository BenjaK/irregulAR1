

#' Simulate from a stationary Gaussian AR(1) process.
#'
#' Simulate from a stationary Gaussian AR(1) process, either at consecutive or
#' irregular times.
#' @param times The time points to simulate for. If a single positive integer,
#'   then that many consecutive values will be simulated.
#' @param rho A real number strictly less than 1 in absolute value.
#' @param sigma A positive real number.
#' @param mu A vector of expected values with length \code{length(times)}, or a
#'   scalar (default equal to 0).
#' @return If \code{length(times) > 1}, a vector of length \code{length(times)}.
#'   Else (\code{length(times) == 1} a vector with as many values as the integer
#'   \code{times}.
#' @export
#' @examples
#' times <- c(3, 5:7, 10)
#' rho <- 0.5
#' sigma <- 1
#' mu <- seq_along(times)
#' ar1_sim(times[1], rho, sigma)
#' ar1_sim(times, rho, sigma)
#' ar1_sim(times, rho, sigma, mu)
ar1_sim <- function(times, rho, sigma, mu = 0) {
  if (length(sigma) != 1 || length(rho) != 1) {
    stop("sigma and rho must be scalars.")
  }
  if (sigma < 0) stop("sigma must be positive")
  if (abs(rho) >= 1) stop("rho must be less than 1 in magnitude.")
  if (length(times) == 1 && times < 0) stop("times must be non-negative.")
  if (length(times) == 1 && times == 0) return(numeric(0))
  if (length(times) == 1 && length(mu) > 1 && length(mu) != times)
    stop("mu should contain as many elements as the number of timepoints.")
  if (length(mu) > 1 && length(times) > 1 && length(mu) != length(times))
    stop("Lengths of times and mu (if not a single value) must match.")
  if (!all(times == as.integer(times))) stop("times must be an integer vector.")

  if (length(times) == 1) {
    return(as.vector(mu + ar1_sim_cpp(times, rho, sigma)))
  } else {
    return(as.vector(mu + ar1_sim_irregular_cpp(times, rho, sigma)))
  }
}

#' Simulate from a stationary Gaussian AR(1) process.
#'
#' Simulate from a stationary Gaussian AR(1) process at \code{n} consecutive
#' time points.
#' @param pred_times A vector of time points to simulate at.
#' @param obs_times A vector of time points at which observations have been
#'   made.
#' @param x_obs The observed values of the process.
#' @param rho A real number strictly less than 1 in absolute value.
#' @param sigma A positive real number.
#' @param mu_pred A vector or scalar with expected values.
#' @param mu_obs A vector or scalar with expected values.
#' @return A vector of length \code{length(pred_times)} with the process
#'   values.
#' @export
#' @examples
#' t_pred <- c(1, 3, 6:8, 10)
#' t_obs <- c(2, 5, 11:12)
#' x_obs <- rnorm(4)
#' rho <- 0.5
#' sigma <- 1
#' # Means equal 0
#' ar1_sim_conditional(t_pred, t_obs, x_obs, rho, sigma)
#' # Time-varying means
#' mu_pred <- t_pred + rnorm(length(t_pred))
#' mu_obs <- t_obs + rnorm(length(t_obs))
#' ar1_sim_conditional(t_pred, t_obs, x_obs + mu_obs, rho, sigma,
#'                     mu_pred, mu_obs)
ar1_sim_conditional <- function(pred_times, obs_times, x_obs, rho, sigma,
                                mu_pred = 0, mu_obs = 0) {
  if (length(sigma) != 1 || length(rho) != 1) {
    stop("sigma and rho must be scalars.")
  }
  if (abs(rho) >= 1) stop("rho must be less than 1 in magnitude.")
  if (sigma < 0) stop("sigma must be positive")
  if (!(length(mu_pred) %in% c(1, length(pred_times)))) {
    stop("mu_pred must be of length 1 or length(pred_times)")
  }
  if (!(length(mu_obs) %in% c(1, length(obs_times)))) {
    stop("mu_obs must be of length 1 or length(obs_times)")
  }
  if (length(mu_pred) == 1) {
    mu_pred <- rep(mu_pred, length(pred_times))
  }
  if (length(mu_obs) == 1) {
    mu_obs <- rep(mu_obs, length(obs_times))
  }
  return(as.vector(ar1_sim_conditional_cpp(pred_times, mu_pred,
                                           x_obs, obs_times, mu_obs,
                                           rho, sigma)))
}