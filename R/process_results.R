#' Evaluate a policy
#'
#' This function evaluates the optimal policy derived from \code{theta} and gives the upper bound of the constraint estimator. 
#' It updates \code{mu0} and \code{nu0} following the estimation step from the alternating optimization procedure. This enables 
#' targeted estimation of the objective functions: risk, constraint, and the main objective, providing a consistent upper bound 
#' for the constraint estimator.
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, 
#' used to recover an extremal point for convex function construction.
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#' @param A A binary vector or matrix of length n indicating treatment assignment 
#' (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes 
#' (in `[0,1]`).
#' @param Xi A numeric vector or matrix of length n indicating adverse events 
#' (0 or 1).
#' @param mu0 A fold-specific function predicting primary outcome (Y) given 
#' treatment (A) and covariates (X).
#' @param nu0 A fold-specific function predicting adverse event outcome (Xi) 
#' given treatment (A) and covariates (X).
#' @param prop_score A function that estimates the propensity score given 
#' treatment (A) and covariates (X).
#' @param lambda A non-negative numeric scalar controlling the penalty for 
#' violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in 
#' `[0,1/2]`, 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the 
#' probability function.
#' @param centered A logical value indicating whether to apply centering in 
#' \code{sigma_beta} (FALSE by default).
#' @param conf_level A numeric scalar representing the confidence level (0.05 
#' by default).
#' @param card_lambda A numeric scalar representing the cardinal of Lambdas.
#' @param card_B A numeric scalar representing the cardinal of B
#' @param pi_X A numeric vector of treatment recommendations associated with 
#' \code{X} (length n). Defaults to NULL. 
#'
#' @return A vector of results for given `theta`.
#' @export
process_results <- function(theta, X, A, Y, Xi, mu0, nu0, prop_score, lambda, 
                            alpha=0.1,  beta=0.05, centered=FALSE, 
                            conf_level = 0.05, card_lambda, card_B, pi_X = NULL) {
  # Correct estimators
  offset_mu <- stats::qlogis(mu0(A,X))
  offset_nu <- stats::qlogis(nu0(A,X))
  H_XA <- HX(A, X, prop_score)
  
  if(is.null(pi_X)){
    psi<- make_psi(theta)
    psi_X <- psi(X)
    sigma_psi_X <- sigma_beta(psi_X,beta, centered) 
    
    # Build targeted risk estimator
    df_mu <- tibble::tibble(
      y = Y, new.cov=H_XA*as.vector(psi_X))
    
    mu_update_obj <- stats::glm(y ~ -1 + ., offset=offset_mu, data = df_mu, 
                                family=stats::binomial())
    epsilon1 <- as.matrix(as.numeric(mu_update_obj$coefficients))
    
    Delta_mu <- function(X) { update_mu_XA(stats::qlogis(mu0(rep(1,nrow(X)),X)), 
                                           epsilon1, psi_X, HX(rep(1,nrow(X)),
                                                               X, prop_score)) - 
        update_mu_XA(stats::qlogis(mu0(rep(0,nrow(X)),X)), epsilon1, psi_X, 
                     HX(rep(0,nrow(X)), X, prop_score)) }
    r_value <- R_p(psi, X, Delta_mu)
    
    # Build targeted constraint estimator
    df_nu <- tibble::tibble(
      xi = Xi, new.cov=H_XA*as.vector(sigma_psi_X))
    
    nu_update_obj <- stats::glm(xi ~ -1 + ., offset=offset_nu, data = df_nu, 
                                family=stats::binomial())
    epsilon2 <- as.matrix(as.numeric(nu_update_obj$coefficients))
    
    Delta_nu <- function(X) { update_nu_XA(stats::qlogis(nu0(rep(1,nrow(X)),X)), epsilon2, sigma_psi_X, HX(rep(1,nrow(X)),X,prop_score)) - 
        update_nu_XA(stats::qlogis(nu0(rep(0,nrow(X)),X)), epsilon2, 
                     sigma_psi_X, HX(rep(0,nrow(X)), X, prop_score))}
    
    constraint <- S_p(psi, X, beta, alpha, centered, Delta_nu)
    l_value <- Lagrangian_p(psi, X, Delta_mu, Delta_nu, lambda, alpha, beta, 
                            centered)
    pi_X <- stats::rbinom(nrow(X), 1, prob = sigma_psi_X)
  }else{
    sigma_psi_X <- pi_X
    df_nu <- tibble::tibble(xi = Xi, new.cov = H_XA * as.vector(sigma_psi_X))
    nu_update_obj <- stats::glm(xi ~ -1 + ., offset = offset_nu, 
                                data = df_nu, family = stats::binomial())
    epsilon2 <- as.matrix(as.numeric(nu_update_obj$coefficients))
    
    Delta_nu <- function(X) {
      update_nu_XA(stats::qlogis(nu0(rep(1, nrow(X)), X)), epsilon2, 
                   sigma_psi_X, HX(rep(1, nrow(X)), X, prop_score)) - 
        update_nu_XA(stats::qlogis(nu0(rep(0, nrow(X)), X)), epsilon2, 
                     sigma_psi_X, HX(rep(0, nrow(X)), X, prop_score))
    }
    r_value <- NA
    constraint <- mean(pi_X * Delta_nu(X)) - alpha
    l_value <- NA
  }  
  
  quantile_V <- stats::qnorm(1-(conf_level/(card_B*card_lambda)))
  epsilon_model <- stats::glm(Y ~ -1 + offset(mu0(pi_X, X)) + 
                                ifelse(pi_X == A, HX(pi_X, X, prop_score), 0), 
                              family = stats::gaussian())
  V_pn <- mean(mu0(pi_X, X) + stats::coef(epsilon_model) * ifelse(pi_X == A, 
                                                           HX(pi_X, X, 
                                                              prop_score), 0))
  Var_pn <- stats::var((ifelse(A == pi_X, 1, 0)/prop_score(A, X)) * (Y - mu0(A, X) + 
                                                                mu0(pi_X, X)) 
                - V_pn)
  
  lower_bound_V_pn <- V_pn - quantile_V * sqrt(Var_pn/nrow(X))
  
  # If policy is 0 everywhere, no need to target estimator of constraint
  if(all(sigma_psi_X==0)){
    results <- data.frame(
      lambda = lambda,
      beta = beta,
      risk = r_value,
      constraint = -alpha,
      obj = l_value,
      policy_value= V_pn,
      lwr_bound_policy_value = lower_bound_V_pn, 
      upper_bound_constraint = -alpha)
    colnames(results) <- c("lambda","beta","risk","constraint","obj", "policy_value", "lwr_bound_policy_value", "upper_bound_constraint")
    return(results)
    }
  
  quantile_S <- stats::qnorm(1-(conf_level/card_B))
  updated_nuXA <- update_nu_XA(stats::qlogis(nu0(A,X)), epsilon2, 
                               sigma_psi_X, H_XA)
  
  Vs_n <- stats::var(H_XA * sigma_psi_X * (Xi - updated_nuXA) + sigma_psi_X * 
                Delta_nu(X) - constraint)
  upper_bound_sn <- constraint + quantile_S * sqrt(Vs_n/nrow(X))
  
  # Extract the policy for the current index
  results <- data.frame(lambda = lambda, beta = beta, risk = r_value, 
                        constraint = constraint, obj = l_value, 
                        policy_value = V_pn, 
                        lwr_bound_policy_value = lower_bound_V_pn, 
                        upper_bound_constraint = upper_bound_sn)
  colnames(results) <- c("lambda","beta","risk","constraint","obj", "policy_value", "lwr_bound_policy_value", "upper_bound_constraint")
  return(results) # Return the updated results for this index
}


#' Oracular evaluation of a policy
#'
#' This function evaluates the optimal policy derived from \code{theta}. This enables the approximation of the objective 
#' functions: risk, constraint, and the main objective and policy value. 
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, used to recover an extremal point for convex function construction.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in `[0,1/2]`, 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param policy A function of \code{X}. Defaults to NULL
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
oracular_process_results <- function(theta, ncov=10L, lambda, alpha=0.1, 
                                     beta=0.05, centered=FALSE, policy = NULL) {
  exp<- generate_data(n=n,ncov=ncov)
  df_obs <- exp[[2]]
  X <- df_obs |>
    dplyr::select(tidyselect::starts_with("X."))|> as.matrix()
  delta_Mu <- exp[[3]]
  delta_Nu <- exp[[4]]
  
  if(is.null(policy)){
    psi<- make_psi(theta)
    r_value <- R_p(psi=psi, X, delta_Mu)
    constraint_value <- S_p(psi=psi, X=X, beta=beta, alpha=alpha, centered=centered, delta_Nu)
    obj_value <- Lagrangian_p(psi, X, delta_Mu, delta_Nu, lambda, alpha, beta, centered)
    # Compute optimal policy value and its lower bound
    Value_policy <- V_p(psi, beta=beta, centered=centered, alpha=alpha, ncov=ncov)
  }else{
    pi_X <- policy(X)
    r_value <- obj_value <- NA 
    constraint_value <- binary_S_p(pi_X, alpha=alpha, delta_Nu = delta_Nu)
    Value_policy <- V_P(psi=NULL, beta=beta, centered=centered, 
                        alpha=alpha, ncov=ncov, policy = policy)
  }
  
  # Extract the policy for the current index
  results <- data.frame(
    lambda = lambda,
    beta = beta,
    risk = r_value,
    constraint = constraint_value,
    obj = obj_value,
    policy_value= Value_policy)
  colnames(results) <- c("lambda","beta","risk","constraint","obj", "policy_value")
  return(results) # Return the updated results for this index
}

#' Select Optimal Beta and Lambda Combination
#'
#' This function loads intermediate results corresponding to lambda-optimal solutions for each beta value.
#' It identifies and returns the beta-lambda combination that minimizes the objective function.
#'
#' @param combinations A matrix or data frame where each row corresponds to a beta and optimal-lambda pair.
#' @param root.path Path to the folder where all results are to be saved.
#'
#' @return A vector of intermediate results including the optimal: lambda, beta, risk, constraint, obj. 
#' @export
get_opt_beta_lambda <- function(combinations,root.path){
  target_filenames <- paste0(combinations[,1], "_", combinations[,2], ".rds")
  all_files <- list.files(
    file.path(root.path, "Evaluation"), 
    pattern = "\\.rds$", full.names = TRUE)
  matched_files <- all_files[basename(all_files) %in% target_filenames]
  optimal_solutions <- do.call(rbind,lapply(matched_files, readRDS))
  optimal_combination <- optimal_solutions[which.max(optimal_solutions$lwr_bound_policy_value),]
  return(optimal_combination)
}
