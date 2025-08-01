% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim_functions.R
\name{FW}
\alias{FW}
\title{Frank-Wolfe algorithm}
\usage{
FW(
  X,
  delta_Mu,
  delta_Nu,
  lambda,
  alpha = 0.1,
  beta = 0.05,
  centered = FALSE,
  precision = 0.05,
  verbose = TRUE
)
}
\arguments{
\item{X}{A matrix of covariates of size n x d (input data).}

\item{delta_Mu}{A function of \code{X} that determines the contrast between primary outcomes.}

\item{delta_Nu}{A function of \code{X} that determines the contrast between adverse event outcomes.}

\item{lambda}{A non-negative numeric scalar controlling the penalty for violating the constraint.}

\item{alpha}{A numeric scalar representing the constraint tolerance (in \link{0,1/2}, 0.1 by default).}

\item{beta}{A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).}

\item{centered}{A logical (FALSE by default) indicating whether to center the policy.}

\item{precision}{A numeric scalar defining the desired convergence precision (0.05 by default). The number of Frank-Wolfe iterations (K) is inversely proportional to this value, calculated as 1/precision.}

\item{verbose}{A logical value indicating whether to print progress updates. Default is \code{TRUE}.}
}
\value{
A numeric matrix containing the optimized parameter \code{theta},
where each row represents the k-th \code{theta} solution at iteration \code{k}.
}
\description{
Implements the Frank-Wolfe optimization algorithm to iteratively refine a convex
combination function \code{psi}. At each iteration, a new solution \code{theta}
is computed via stochastic gradient descent (SGD) and added to the convex combination
in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
}
