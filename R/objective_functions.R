#' Risk Function for Conditional Average Treatment Effect (CATE)
#'
#' Computes the risk function \eqn{R_p} for estimating the Conditional Average Treatment Effect (CATE).  
#' The function minimizes the squared error between \code{psi(X)} and \code{delta_Y(X)}.
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix of covariates of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the contrast between primary outcomes.
#'
#' @return A numeric scalar representing the risk function value.
#' @export
R_p <- function(psi, X, delta_Mu){
    out <- mean(psi(X)^2 - 2 * psi(X)* delta_Mu(X))
    return(out)
}

#' Constraint Function
#'
#' Computes the constraint function \eqn{S_p}, which ensures that the learned policy satisfies a constraint. 
#' This function enforces a limit on the expected impact of treatment via \code{delta_Z}.
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix of covariates of size n x d (input data).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param delta_Nu A function of \code{X} that determines the contrast between adverse event outcomes.
#'
#' @return A numeric scalar representing the constraint function value.
#' @export
S_p <- function(psi, X, beta, alpha, centered, delta_Nu){
  psi_X <- psi(X)
  out <- mean(sigma_beta(psi_X, beta, centered) * delta_Nu(X)) - alpha
    return(out)
}

#' Objective Function taking the form of a Lagrangian
#'
#' Computes the objective function, which balances the risk function \code{R_p}  
#' and the constraint function \code{S_p} using a parameter \code{lambda}.
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix of covariates of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the contrast between primary outcomes.
#' @param delta_Nu A function of \code{X} that determines the contrast between adverse event outcomes.
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A numeric scalar representing the objective function value.
#' @export
Lagrangian_p <- function(psi, X, delta_Mu, delta_Nu, lambda, alpha=0.1, beta=0.05, centered=FALSE){
    out <- R_p(psi, X,delta_Mu) + lambda*S_p(psi, X, beta, alpha, centered, delta_Nu)
    return(out)
}

#' Gradient of the Objective Function
#'
#' Computes the gradient of the objective function with respect to \code{psi} at X.  
#' The gradient is used in optimization algorithms like Stochastic Gradient Descent (SGD).
#'
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix of covariates of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the contrast between primary outcomes.
#' @param delta_Nu A function of \code{X} that determines the contrast between adverse event outcomes.
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A numeric vector representing the gradient of the objective function with respect to \code{psi(X)}.
#' @export
grad_Lagrangian_p <- function(psi, X, delta_Mu, delta_Nu, lambda, alpha=0.1, beta=0.05, centered=FALSE){
  psi_X <- psi(X)
  2*(psi_X - delta_Mu(X)) + lambda*sigma_beta_prime(psi_X, beta, centered)*delta_Nu(X)
}

#' @rdname grad_Lagrangian_p
#' 
#' @export
grad_Lagrangian_p_X <- function(psi_X, delta_Mu_X, delta_Nu_X, lambda, alpha=0.1, beta=0.05, centered=FALSE){
  2*(psi_X - delta_Mu_X) + lambda*sigma_beta_prime(psi_X, beta, centered)*delta_Nu_X
}

