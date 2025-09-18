#' Learn Optimal Decision Threshold
#'
#' This function estimates the optimal decision threshold for treatment assignment 
#' using cross-fitted nuisance function estimates and targeted maximum likelihood 
#' estimation (TMLE). The procedure evaluates candidate thresholds based on empirical 
#' performance criteria, selecting the threshold that maximizes the policy value 
#' under constraint satisfaction.
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, used to recover an extremal point for convex function construction.
#' @param X A matrix or data frame of covariates of size n x d (input data in [0,1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of adverse events outcomes.
#' @param folds A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}. 
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#'
#' @return 
#' A numeric value corresponding to the optimal threshold chosen 
#' from the candidate sequence.
#' 
#' @export
learn_threshold <- function(theta, X, A, Y, Xi, folds, alpha){
  check_data(Y, Xi, A, X, folds)
  
  X_train <- X[folds[[2]],]
  A_train <- A[folds[[2]]]
  Y_train <- Y[folds[[2]]]
  Xi_train <- Xi[folds[[2]]]

  mu.hat.nj <- PLUCR::estimate_mu(Y=Y, A=A, X=X, folds=folds, SL.library = SL.library, V=2L)
  nu.hat.nj <- PLUCR::estimate_nu(Xi=Xi, A=A, X=X, folds=folds, SL.library = SL.library,V=2L)
  ps.hat.nj <- PLUCR::estimate_ps(A=A, X=X, folds=folds, SL.library = SL.library,V=2L)
  
  mu0_train <- function(a,x){mu.hat.nj(a=a,x=x,ff=1)}
  nu0_train <- function(a,x){nu.hat.nj(a=a,x=x,ff=1)}
  prop_score_train <- function(a,x){ps.hat.nj(a,x,1)}
  
  sigma_psi <- sigma_beta(make_psi(theta)(X), beta=attr(theta, "beta"))
  
  thresholds <- seq(0,1, 5*1e-2)
  results <- NULL
  for(t in thresholds){
    sample_id <- sample(c(folds[[1]],folds[[2]]),size = 2*length(folds[[2]])/3)
    X_t <- X[sample_id,]
    A_t <- A[sample_id]
    Y_t <- Y[sample_id]
    Xi_t <- Xi[sample_id]
    pi_train <- ifelse(sigma_psi[sample_id]>t, 1,0)
    corr_nparams <- lwr_upper_bound_estimators(mu0=mu0_train, nu0=nu0_train, prop_score=prop_score_train, pi=pi_train, X=X_t, A=A_t, Y=Y_t, Xi=Xi_t, alpha=alpha)[[1]]
    results <- rbind(results, c(corr_nparams[[1]], corr_nparams[[2]]))
  }
  final_threshold <- thresholds[which(results[, 2] < 0 & results[, 1] == max(results[results[, 2] < 0, 1]))]
  return(final_threshold)
}


#' Lower and upper bound estimators for policy value and constraint
#'
#' This function computes lower and upper confidence bounds for the policy value 
#' and constraint, respectively, based on targeted maximum likelihood estimation 
#' (TMLE) updates. 
#'
#' @param mu0 A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X).
#' @param nu0 A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X).
#' @param prop_score A function that estimates the propensity score given treatment (A) and covariates (X).
#' @param pi A binary treatment rule vector indicating treatment assignment 
#'   under the learned decision rule.
#' @param X A matrix or data frame of covariates of size n x d (input data in [0,1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of adverse events outcomes.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#'
#' @return 
#' A list containing a numeric vector of length 2:
#' \itemize{
#'   \item \code{[1]}: Lower bound for the primary outcome effect.
#'   \item \code{[2]}: Upper bound for the adverse event effect.
#' }
#' 
#' @export
lwr_upper_bound_estimators <- function(mu0, nu0, prop_score, pi, X, A, Y, Xi, alpha){
  offset_nu <- qlogis(nu0(A, X))
  H_XA <- HX(A, X, prop_score)
  
  df_nu <- tibble::tibble(xi = Xi, new.cov = ifelse(H_XA * as.vector(pi) == 0, 1e-5, H_XA * as.vector(pi)))
  nu_update_obj <- stats::glm(xi ~ -1 + ., offset = offset_nu, 
                              data = df_nu, family = binomial())
  epsilon2 <- as.matrix(as.numeric(nu_update_obj$coefficients))
  Delta_nu <- function(X) {
    nu1_updated <- plogis(qlogis(nu0(rep(1,nrow(X)), X)) + as.vector(as.numeric(epsilon2) *pi) * HX(rep(1,nrow(X)), X, prop_score))
    nu0_updated <- plogis(qlogis(nu0(rep(0,nrow(X)), X)) + as.vector(as.numeric(epsilon2) *pi) * HX(rep(0,nrow(X)), X, prop_score))
    return(nu1_updated - nu0_updated)
  }
  updated_nuXA <- update_nu_XA(qlogis(nu0(A, X)),epsilon2, sigma_psi_collection=pi, H_XA = HX(A, X, prop_score))
  s_n <- binary_S_p(pi, X, alpha, Delta_nu)
  
  Vs_n <- var(ifelse(A==pi,H_XA,0) * (Xi - updated_nuXA) + pi * Delta_nu(X) - s_n)
  upper_bound_sn <- s_n + 1.64 * sqrt(Vs_n/nrow(X))
  
  epsilon_model <- stats::glm(Y ~ -1 + offset(mu0(pi, X)) + ifelse(pi == A, HX(pi, X, prop_score),0), family = gaussian())
  V_targeted<- mean(mu0(pi,X) + coef(epsilon_model) * ifelse(pi==A,HX(pi,X,prop_score),0))
  Var_pn <- var( (ifelse(A==pi,1,0)/prop_score(A,X))*(Y- mu0(A,X) + mu0(pi,X))-V_targeted) #varphi
  lower_bound_V_pn <- V_targeted - 1.64*sqrt(Var_pn/nrow(X))
  return(list(c(lower_bound_V_pn, upper_bound_sn)))
}

#' Naive lower and upper bound estimators for policy value and constraint
#'
#' This function computes lower and upper confidence bounds for the policy value 
#' and constraint, respectively, based on their plug-in estimators.
#'
#' @param mu0 A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X).
#' @param nu0 A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X).
#' @param prop_score A function that estimates the propensity score given treatment (A) and covariates (X).
#' @param pi A numeric vector of binary decisions of length n.
#' @param X A matrix or data frame of covariates of size n x d (input data in [0,1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of adverse events outcomes.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#'
#' @return 
#' A list containing a numeric vector of length 2:
#' \itemize{
#'   \item \code{[1]}: Lower bound for the primary outcome effect.
#'   \item \code{[2]}: Upper bound for the adverse event effect.
#' }
#' 
#' @export
naive_lwr_upper_bound_estimators <- function(mu0, nu0, pi, X, A, Y, Xi, prop_score, alpha){
  V_pn <- V_Pn(pi, mu0(rep(1,nrow(X)),X), mu0(rep(0,nrow(X)),X))
  Var_pn <- var((ifelse(A==pi,1,0)/prop_score(A,X))*(Y- mu0(A,X) + mu0(pi,X))-V_pn) #varphi
  lower_bound_V_pn <- V_pn - 1.64*sqrt(Var_pn/nrow(X))
  Delta_nu <- function(X){
    nu0(rep(1,nrow(X)),X) - nu0(rep(0,nrow(X)),X)
  }
  s_n <- binary_S_p(pi, X, alpha, Delta_nu)
  Vs_n <- var( ifelse(A==pi, HX(A, X, prop_score), 0 ) * (Xi - ifelse(A==1, nu0(rep(1, nrow(X)),X), nu0(rep(0, nrow(X)),X))) + pi*Delta_nu(X)  - s_n)
  upper_bound_sn <- s_n + 1.64 * sqrt(Vs_n/nrow(X))
  return(list(c(lower_bound_V_pn, upper_bound_sn)))
}