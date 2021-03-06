% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dcov_drho}
\alias{dcov_drho}
\title{Derivative of the covariance matrix for a stationary Gaussian AR(1) process.}
\usage{
dcov_drho(times, rho, sigma)
}
\arguments{
\item{times}{An vector of positive integers, preferably ordered.}

\item{rho}{A real number strictly less than 1 in absolute value.}

\item{sigma}{A positive real number.}
}
\value{
A square matrix with \code{length(times)} rows.
}
\description{
Creates the derivate of the covariance matrix of an AR(1) process with
respect to the parameter \code{rho}. The process has been observed at the
time points in the vector \code{times} and is assumed to be in stationarity,
and to have Gaussian errors.
}
\examples{
times <- c(1, 4:5, 7)
rho <- 0.5
sigma <- 1
dcov_drho(times, rho, sigma)
}
