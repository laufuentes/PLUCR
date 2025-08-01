% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sigma_beta_prime}
\alias{sigma_beta_prime}
\title{Derivative of link function}
\usage{
sigma_beta_prime(t, beta = 0.05, centered = FALSE)
}
\arguments{
\item{t}{A vector of numerics (in \link{-1,1}).}

\item{beta}{A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).}

\item{centered}{A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).}
}
\value{
The derivative of \code{sigma_beta} evaluated at t.
}
\description{
Computes the derivative of the link function \code{sigma_beta},
with respect to t.
}
