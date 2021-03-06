% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ar1_cov_consecutive}
\alias{ar1_cov_consecutive}
\title{Covariance matrix for a stationary Gaussian AR(1) process, observed at
consecutive timepoints.}
\usage{
ar1_cov_consecutive(n, rho, sigma)
}
\arguments{
\item{n}{An integer greater than or equal to 1.}

\item{rho}{A real number strictly less than 1 in absolute value.}

\item{sigma}{A positive real number.}
}
\value{
A matrix with \code{n} rows and \code{n} columns.
}
\description{
Creates the covariance matrix of an AR(1) process with parameters \code{rho}
and \code{sigma}, observed at \code{n} consecutive time points. The process
is assumed to be in stationarity and to have Gaussian errors.
}
\keyword{internal}
