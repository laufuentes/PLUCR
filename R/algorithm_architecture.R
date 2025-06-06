#' Evaluate a policy
#'
#' This function estimates a final policy by aggregating fold-specific contrast
#' functions (`Delta_mu` and `Delta_nu`) using cross-validation. The contrasts
#' are applied fold-wise, and the final policy parameters are estimated using
#' a policy learning algorithm (e.g., `FW`).
#'
#' @param theta Covariate matrix (n x p) for all units.
#' @param X Covariate matrix (n x p) for all units.
#' @param delta_Mu A list of fold-specific treatment contrast functions (e.g., T-learner estimates).
#' @param delta_Nu A list of fold-specific constraint contrast functions.
#' @param lambda Regularization parameter for the risk component.
#' @param alpha Constraint tolerance (typically between 0 and 1).
#' @param beta Regularization parameter for the risk component.
#' @param centered Logical; whether to center the features.
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
process_results <- function(theta, X, A, Y, Xi, mu0, nu0, prop_score, lambda, alpha,  beta, centered) {
  # Correct estimators
  offset_mu <- qlogis(mu0(A,X))
  offset_nu <- qlogis(nu0(A,X))
  
  psi<- make_psi(theta)
  psi_X <- psi(X)
  sigma_psi_X <- sigma_beta(psi_X,beta, centered)
  H_XA <- HX(A, X, prop_score)
  
  df_mu <- tibble::tibble(
    y = Y, new.cov=H_XA*as.vector(psi_X))
  
  df_nu <- tibble::tibble(
    xi = Xi, new.cov=H_XA*as.vector(sigma_psi_X))
  
  mu_update_obj <- stats::glm(y ~ -1 + ., offset=offset_mu, data = df_mu, family=binomial())
  epsilon1 <- as.matrix(as.numeric(mu_update_obj$coefficients))
  
  nu_update_obj <- stats::glm(xi ~ -1 + ., offset=offset_nu, data = df_nu, family=binomial())
  epsilon2 <- as.matrix(as.numeric(nu_update_obj$coefficients))
  
  Delta_mu <- function(X) { update_mu_XA(qlogis(mu0(rep(1,nrow(X)),X)), epsilon1, psi_X, HX(rep(1,nrow(X)),X,prop_score)) - 
      update_mu_XA(qlogis(mu0(rep(0,nrow(X)),X)), epsilon1, psi_X, HX(rep(0,nrow(X)),X,prop_score)) }
  Delta_nu <- function(X) { update_nu_XA(qlogis(nu0(rep(1,nrow(X)),X)), epsilon2, sigma_psi_X, HX(rep(1,nrow(X)),X,prop_score)) - 
      update_nu_XA(qlogis(nu0(rep(0,nrow(X)),X)), epsilon2, sigma_psi_X, HX(rep(0,nrow(X)),X,prop_score)) }
  # Extract the policy for the current index
  results <- data.frame(
    lambda = lambda,
    beta = beta,
    risk = R_p(psi, X, Delta_mu),
    constraint = S_p(
      psi,X,
      beta,alpha, centered, 
      Delta_nu),
    obj = Lagrangian_p(psi, X, Delta_mu, Delta_nu, lambda, alpha, beta, centered))
  colnames(results) <- c(
    "lambda",
    "beta",
    "risk",
    "constraint",
    "obj")
  updated_nuXA <- update_nu_XA(qlogis(nu0(A,X)), epsilon2, sigma_psi_X, H_XA)
  V_n <- var(H_XA* sigma_psi_X*(Xi- updated_nuXA) +
               sigma_psi_X* Delta_nu(X) -results$constraint)
  upper_bound <- results$constraint + 1.64*sqrt(V_n/nrow(X))
  
  return(list(results, upper_bound)) # Return the updated results for this index
}

#' Learn Final Policy Using Cross-Validated Contrast Estimators
#'
#' This function estimates a final policy by aggregating fold-specific contrast
#' functions (`Delta_mu` and `Delta_nu`) using cross-validation. The contrasts
#' are applied fold-wise, and the final policy parameters are estimated using
#' a policy learning algorithm (e.g., `FW`).
#'
#' @param lambda Regularization parameter for the risk component.
#' @param beta Regularization parameter for the constraint component.
#' @param X Covariate matrix (n x p) for all units.
#' @param s A vector indicating fold assignments (length n).
#' @param Delta_mu_nj_folds A list of fold-specific treatment contrast functions (e.g., T-learner estimates).
#' @param Delta_nu_nj_folds A list of fold-specific constraint contrast functions.
#' @param alpha Constraint tolerance (typically between 0 and 1).
#' @param centered Logical; whether to center the features.
#' @param precision A numeric scalar that determines the convergence precision desired.
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
Final_policy <- function(lambda, beta,X, s, Delta_mu_nj_folds, Delta_nu_nj_folds, alpha, centered, precision){
  Delta_mu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      delta_mu <- Delta_mu_nj_folds[[fold]]
      out[s==fold] = delta_mu(X)
    }
    return(out)
  }
  Delta_nu_CV <- function(X){
    out <- rep(0,nrow(X))
    for(fold in unique(s)){
      X_fold <- X[s==fold,]
      delta_nu <- Delta_nu_nj_folds[[fold]]
      out[s==fold] = delta_nu(X)
    }
    return(out)
  }
  theta_final <- FW(X, Delta_mu_CV, Delta_nu_CV, lambda, alpha, 
                    beta, centered, precision)
  return(theta_final)
}