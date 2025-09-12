#' Evaluate a policy
#'
#' This function evaluates the optimal policy derived from \code{theta} and gives the upper bound of the constraint estimator. 
#' It updates \code{mu0} and \code{nu0} following the estimation step from the alternating optimization procedure. This enables 
#' targeted estimation of the objective functions: risk, constraint, and the main objective, providing a consistent upper bound 
#' for the constraint estimator.
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, used to recover an extremal point for convex function construction.
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of length n indicating adverse events (0 or 1).
#' @param mu0 A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X).
#' @param nu0 A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X).
#' @param prop_score A function that estimates the propensity score given treatment (A) and covariates (X).
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
process_results <- function(theta, X, A, Y, Xi, mu0, nu0, prop_score, lambda, alpha=0.1,  beta=0.05, centered=FALSE) {
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
      update_nu_XA(qlogis(nu0(rep(0,nrow(X)),X)), epsilon2, sigma_psi_X, HX(rep(0,nrow(X)),X,prop_score))}
  
  # Compute optimal policy value and its lower bound
  pi_X <- rbinom(nrow(X),1, prob= sigma_psi_X) # binary policy
  epsilon_model <- stats::glm(Y ~ -1 + offset(mu0(pi_X, X)) + ifelse(pi_X == A, HX(pi_X, X, prop_score), 0),
    family = gaussian())
  
  V_pn <- mean(mu0(pi_X, X) + 
                 coef(epsilon_model) * ifelse(pi_X == A, HX(pi_X, X, prop_score), 0))
  
  
  Var_pn <- var( (ifelse(A==pi_X,1,0)/prop_score(A,X))*(Y- mu0(A,X) + mu0(pi_X,X))-V_pn) #varphi
  upper_bound_V_pn <- V_pn - 1.64*sqrt(Var_pn/nrow(X))
  
  # Extract the policy for the current index
  results <- data.frame(
    lambda = lambda,
    beta = beta,
    risk = R_p(psi, X, Delta_mu),
    constraint = S_p(
      psi, X,
      beta, alpha, centered, 
      Delta_nu),
    obj = Lagrangian_p(psi, X, Delta_mu, Delta_nu, lambda, alpha, beta, centered),
    policy_value= V_pn,
    lwr_bound_policy_value = upper_bound_V_pn)
  colnames(results) <- c("lambda","beta","risk","constraint","obj", "policy_value", "lwr_bound_policy_value")
  
  # Compute upper bound for constraint
  updated_nuXA <- update_nu_XA(qlogis(nu0(A,X)), epsilon2, sigma_psi_X, H_XA)
  
  Vs_n <- var(H_XA* sigma_psi_X*(Xi- updated_nuXA) +
                sigma_psi_X* Delta_nu(X) -results$constraint)
  upper_bound_sn <- results$constraint + 1.64*sqrt(Vs_n/nrow(X))
  
  return(list(results, upper_bound_sn)) # Return the updated results for this index
}


#' Naive policy evaluation
#'
#' This function evaluates the optimal policy derived from \code{theta} and gives the upper bound of the constraint estimator using the naive approach, 
#' that is to say the baseline nuisances. 
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, used to recover an extremal point for convex function construction.
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of length n indicating adverse events (0 or 1).
#' @param mu0 A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X).
#' @param nu0 A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X).
#' @param prop_score A function that estimates the propensity score given treatment (A) and covariates (X).
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
naive_process_results <- function(theta, X, A, Y, Xi, mu0, nu0, prop_score, lambda, alpha=0.1,  beta=0.05, centered=FALSE) {
  psi<- make_psi(theta)
  psi_X <- psi(X)
  sigma_psi_X <- sigma_beta(psi_X,beta, centered)
  H_XA <- HX(A, X, prop_score)
  
  Delta_mu <- function(X) { mu0(rep(1,nrow(X)),X)- mu0(rep(0,nrow(X)),X)}
  Delta_nu <- function(X) { nu0(rep(1,nrow(X)),X) - nu0(rep(0,nrow(X)),X)}
  
  mu1_X <- mu0(rep(1,nrow(X)),X)
  mu0_X <- mu0(rep(0,nrow(X)),X)
  
  # Compute optimal policy value and its lower bound
  pi_X <- rbinom(nrow(X),1, prob= sigma_psi_X) # binary policy
  V_pn <- V_Pn(pi_X, mu1_X, mu0_X) # estimated policy value 
  
  mu_AX <- mu0(A,X)
  mu_pi_X <- mu0(pi_X,X)
  
  Var_pn <- var( (ifelse(A==pi_X,1,0)/prop_score(A,X))*(Y- mu_AX + mu_pi_X)-V_pn) #varphi
  lower_bound_V_pn <- V_pn - 1.64*sqrt(Var_pn/nrow(X))
  
  # Extract the policy for the current index
  results <- data.frame(
    lambda = lambda,
    beta = beta,
    risk = R_p(psi, X, Delta_mu),
    constraint = S_p(psi, X,beta, alpha, centered, Delta_nu),
    obj = Lagrangian_p(psi, X, Delta_mu, Delta_nu, lambda, alpha, beta, centered),
    policy_value= V_pn,
    lwr_bound_policy_value = lower_bound_V_pn)
  colnames(results) <- c("lambda","beta","risk","constraint","obj", "policy_value", "lwr_bound_policy_value")
  
  # Compute upper bound for constraint
  nuXA <- nu0(A,X)
  
  Vs_n <- var(H_XA* sigma_psi_X*(Xi- nuXA) +
                sigma_psi_X* Delta_nu(X) -results$constraint)
  upper_bound_sn <- results$constraint + 1.64*sqrt(Vs_n/nrow(X))
  
  return(list(results, upper_bound_sn)) # Return the updated results for this index
}

#' Oracular evaluation of a policy
#'
#' This function evaluates the optimal policy derived from \code{theta}. This enables the approximation of the objective 
#' functions: risk, constraint, and the main objective and policy value. 
#'
#' @param theta A numeric matrix (k x d). Each row is from FW inner minimization, used to recover an extremal point for convex function construction.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param scenario_mu String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix").
#' @param scenario_nu String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix").
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#'
#' @return A vector of optimized policy parameters (`theta`) trained across folds.
#' @export
oracular_process_results <- function(theta, ncov=10L, 
                                     scenario_mu=c("Linear", "Threshold", "Mix", "Null", "Constant", "Realistic"), 
                                     scenario_nu=c("Linear", "Threshold", "Mix", "Satisfied", "Realistic"), 
                                     lambda, alpha=0.1,  beta=0.05, centered=FALSE) {
  psi<- make_psi(theta)
  if(scenario_mu=="Realistic"){
    exp<- generate_realistic_data(1e6)
    df_complete <- exp[[1]]
    df_obs <- exp[[2]]
    X <- df_obs %>%
      dplyr::select(starts_with("X."))%>% as.matrix()
    X_norm <- phi(X)
    attr(X_norm, "min_Y") <- unique(df_complete$min_Y)
    attr(X_norm, "max_Y") <- unique(df_complete$max_Y)
    delta_Mu <- function(x)exp[[3]](phi_inv(x))
    delta_Nu <- function(x)exp[[4]](phi_inv(x))
    X <- X_norm
  }else{
    exp<- generate_data(n=n,ncov=ncov, scenario_mu=scenario_mu, scenario_nu=scenario_nu)
    df_obs <- exp[[2]]
    X <- df_obs %>%
      dplyr::select(starts_with("X."))%>% as.matrix()
    delta_Mu <- exp[[3]]
    delta_Nu <- exp[[4]]
  }
  # Compute optimal policy value and its lower bound
  Value_policy <- V_p(psi, beta=beta, centered=centered, alpha=alpha, ncov=ncov, 
                      scenario_mu=scenario_mu, scenario_nu=scenario_nu)
  
  # Extract the policy for the current index
  results <- data.frame(
    lambda = lambda,
    beta = beta,
    risk = R_p(psi=psi, X, delta_Mu),
    constraint = S_p(psi=psi, X=X, beta=beta, alpha=alpha, centered=centered, delta_Nu),
    obj = Lagrangian_p(psi, X, delta_Mu, delta_Nu, lambda, alpha, beta, centered),
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
