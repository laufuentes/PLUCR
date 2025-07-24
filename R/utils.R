#' The expit/logit functions
#'
#' @description
#' expit(x) returns \code{1/(1+exp(-x))} for any real x. 
#' 
#' logit(p) return \code{log(p/(1-p))} for any probability p.
#' 
#' @export
expit <- plogis

#' @rdname expit
#' @export
logit <- qlogis

#' Check Input Data for Validity
#'
#' Performs quality control checks on the input data to ensure it meets expected formats and conditions.
#'
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of length n indicating adverse events (0 or 1).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param folds A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}. 
#'
#' @return A list with elements: `ok` (logical), and `diagnoses` (character vector of issues).
#' @export
check_data <- function(Y, Xi, A, X, folds) {
  diagnoses <- c()
  ok <- TRUE
  
  # Check Y values
  if (!is.numeric(Y) || !all(Y > 0 & Y < 1)) {
    diagnoses <- c(diagnoses, "Y must be numeric and contain values strictly between 0 and 1.")
    warning("Y must be numeric and contain values strictly between 0 and 1.")
    ok <- FALSE
  }
  
  # Check Xi binary
  if (!all(Xi %in% c(0, 1))) {
    diagnoses <- c(diagnoses, "Xi must be binary (only 0 and 1 values allowed).")
    warning("Xi must be binary (only 0 and 1 values allowed).")
    ok <- FALSE
  }
  
  # Check A binary
  if (!all(A %in% c(0, 1))) {
    diagnoses <- c(diagnoses, "A must be binary (only 0 and 1 values allowed).")
    warning("A must be binary (only 0 and 1 values allowed).")
    ok <- FALSE
  }
  
  # Check X type
  if (!is.matrix(X) && !is.data.frame(X)) {
    # Verifier qu'il n'y a pas des noms interdits
    diagnoses <- c(diagnoses, "X must be a matrix or data frame.")
    warning("X must be a matrix or data frame.")
    ok <- FALSE
  }
  
  if (nrow(X) != length(Y) || nrow(X) != length(Xi) || nrow(X) != length(A)) {
    warning(sprintf("X, Y, Xi and A must have the same number of rows. Got: nrow(X) = %d, length(Y) = %d, length(Xi) = %d, length(A) = %d",
                    nrow(X), length(Y), length(Xi), length(A)))
  }
  
  # Check fold balance
  fold_ids <- length(folds)
  
  for (f in fold_ids) {
    fold_A <- A[folds[[f]]]
    if (!(any(fold_A == 0) && any(fold_A == 1))) {
      msg <- paste("Fold", f, "does not contain both treatment groups (A=0 and A=1).")
      diagnoses <- c(diagnoses, msg)
      warning(msg)
      ok <- FALSE
    }
  }
  
  return(list(ok = ok, diagnoses = diagnoses))
}


#' Link Function
#'
#' Link function mapping \eqn{[-1,1]} to \eqn{[0,1]}, parametrized    
#' by \code{beta} with an optional centering.
#'
#' @param t A vector of numerics (in [-1,1]).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A numeric vector of treatment probabilities.
#' @export
sigma_beta <- function(t, beta=0.05, centered=FALSE) {
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
  if (centered) {
    cent <- 0.5 - c_beta * log(2 / (1 + exp(-beta)))
  } else {
    cent <- 0
  }
  out <- c_beta * log((1 + exp(beta * t)) / (1 + exp(-beta))) + cent
  return(out)
}

#' Derivative of Link Function
#'
#' Computes the derivative of the link function \code{sigma_beta},  
#' with respect to t.
#'
#' @param t A vector of numerics (in [-1,1]).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return The derivative of \code{sigma_beta} evaluated at t.
#' @export
sigma_beta_prime <- function(t, beta=0.05, centered=FALSE){
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  c_beta <- 1 / log(
    (1 + exp(beta)) / (1 + exp(-beta))
  )
  out <- c_beta *(beta*exp(beta*t))/(1+ exp(beta*t))
  return(out)
}


#' Oracular Approximation of Value Function
#'
#' Computes the expected outcome under a policy determined by the previously optimized \code{psi(X)}.  
#' The policy assigns treatment probabilistically based on \code{sigma_beta(psi(X))},  
#' and the expected outcome is calculated using counterfactual outcomes.
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param B Integer, number of Monte Carlo repetitions (1e4 by default).
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param scenario_mu String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix").
#' @param scenario_nu String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix").
#' @param seed Integer or NA (NA by default).
#'
#' @return A numeric scalar representing the expected primary outcome under the policy.
#' @export
V_p <- function(psi, beta=0.05, centered=FALSE, alpha=0.1, B=1e4, ncov=10L, 
                scenario_mu=c("Linear", "Threshold", "Mix"), 
                scenario_nu=c("Linear", "Threshold", "Mix"), seed=NA){
  `%>%`<- magrittr::`%>%`
  df <- generate_data(B, ncov=ncov, scenario_mu=scenario_mu, scenario_nu=scenario_nu, seed=seed)[[1]]
  X <- df%>%dplyr::select(dplyr::starts_with("X."))%>% as.matrix()
  y1 <- df$Y.1
  y0 <- df$Y.0
  
  psi_X <- psi(X)
  sigma_psi <-sigma_beta(psi_X, beta, centered)
  if(!is.na(seed)){
    set.seed(seed)
  }
  action <- stats::rbinom(nrow(X), 1, sigma_psi)
  out <- mean(action * y1 + (1 - action) * y0)
  return(out)
}

#' Estimation of Policy Value
#'
#' Computes the expected outcome under a policy determined by the previously optimized \code{psi(X)}.  
#' The policy assigns treatment probabilistically based on \code{sigma_beta(psi(X))},  
#' and the expected outcome is calculated using counterfactual outcomes.
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param y1 A numeric vector or matrix of length n representing primary outcomes under treatment (in [0, 1]).
#' @param y0 A numeric vector or matrix of length n representing primary outcomes under no treatment (in [0, 1]).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A numeric scalar representing the expected primary outcome under the policy.
#' @export
V_Pn <- function(policy, y1, y0){
  `%>%`<- magrittr::`%>%`
  out <- mean(policy * y1 + (1 - policy) * y0)
  return(out)
}


#' Compute the Inverse Propensity Score Weight
#'
#' This function computes the inverse propensity score weight based on treatment assignment and a propensity score model.
#'
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param prop_score A function that estimates the propensity score given treatment (A) and covariates (X).
#'
#' @return A numeric value representing the inverse propensity score weight.
#'
#' @examples
#' # Example usage:
#' prop_model <- function(A, X) { 0.5 }  # Constant propensity score for illustration
#' HX(1, data.frame(x1 = 1, x2 = 2), prop_model)
#'
#' @export
HX <- function(A, X, prop_score){
  out <- (2*A-1)/prop_score(A,X)
  return(out)
}


