% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_approaches.R
\name{oracular_approach_algorithm}
\alias{oracular_approach_algorithm}
\title{Oracular approach main algorithm}
\usage{
oracular_approach_algorithm(
  X,
  A,
  Y,
  Xi,
  folds,
  delta_Mu,
  delta_Nu,
  scenario_mu,
  scenario_nu,
  Lambdas = seq(1, 8, by = 1),
  alpha = 0.1,
  precision = 0.05,
  B = c(0.05, 0.1, 0.25, 0.5),
  centered = FALSE,
  root.path
)
}
\arguments{
\item{X}{A matrix of covariates of size n x d (input data).}

\item{A}{A binary vector of size n indicating treatment assignment (0 or 1).}

\item{Y}{A numeric vector or matrix of length n representing primary outcomes (in \link{0, 1}).}

\item{Xi}{A numeric vector or matrix of length n indicating adverse events (0 or 1).}

\item{folds}{A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}.}

\item{delta_Mu}{A function that computes the treatment effect (mu difference) from covariates.}

\item{delta_Nu}{A function that computes the selection effect (nu difference) from covariates.}

\item{scenario_mu}{String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix").}

\item{scenario_nu}{String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix").}

\item{Lambdas}{A sequence of non-negative numeric scalars controlling the penalty for violating the constraint (seq(1,8,by=1) by default).}

\item{alpha}{A numeric scalar representing the constraint tolerance (in \link{0,1/2}, 0.1 by default).}

\item{precision}{A numeric scalar defining the desired convergence precision (0.05 by default). The number of Frank-Wolfe iterations (K) is inversely proportional to this value, calculated as 1/precision.}

\item{B}{A vector of non-negative scalars controlling the sharpness of the treatment probability function (c(0.05, 0.1, 0.25, 0.5) by default).}

\item{centered}{A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).}

\item{root.path}{Path to the folder where all results are to be saved.}
}
\value{
A list of matrices (\code{theta_0} and \code{theta_final}), or \code{theta_0} alone. These matrices are used to construct the optimal treatment rule in two steps. First, build \code{psi} using the \code{make_psi} function and evaluate it at \code{X} (i.e., \code{psi(X)}). Then, obtain the optimal treatment rule by applying \code{sigma_beta} to the selected \code{beta} attribute (\code{sigma_beta(psi(X), beta)}).
}
\description{
Learns an optimal treatment policy under constraints using true data structures.
It begins by validating and preprocessing the input data.
Subsequently, it approximates the optimal treatment policy via the Frank-Wolfe optimization algorithm.
The procedure includes an inner grid search over candidate values of \code{lambda} and \code{beta} to identify the policy
that maximizes the expected primary outcome (policy value) while satisfying a constraint on the expected rate of adverse events.
}
