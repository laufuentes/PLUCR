#' Update Mu via Augmented Covariate Adjustment for fixed X
#'
#' Computes updated Mu predictions by adding the bias correction for fixed 
#' \code{psi_collection} scaled by coefficients \code{epsilon1_collection},
#' then transforms back via the logistic function.
#'
#' @param offset_mu_XA A numeric vector of logit-scale baseline mu0 predictions.
#' @param epsilon1 A numeric vector of GLM coefficients for each column in \code{psi_collection}.
#' @param psi_collection A matrix whose columns are optimal \code{psi} solutions.
#' @param H_XA A numeric vector of inverse-propensity weights, typically from \code{HX()}.
#'
#' @return A numeric vector of updated mu on the \[0,1\] scale.
#'
#'@export
update_mu_XA <- function(offset_mu_XA, epsilon1, psi_collection, H_XA){ 
  out <- expit(offset_mu_XA + H_XA*(psi_collection%*%epsilon1)) 
  return(out)
}

#' Update Nu via Augmented Covariate Adjustment for fixed X
#'
#' Computes updated Nu predictions by adding the bias correction for fixed 
#' \code{sigma_psi_collection} scaled by coefficients \code{epsilon2_collection},
#' then transforms back via the logistic function.
#'
#' @param offset_nu_XA A numeric vector of logit-scale baseline nu0 predictions.
#' @param epsilon2 A numeric vector of GLM coefficients for each column in \code{sigma_psi_collection}.
#' @param sigma_psi_collection A matrix whose columns are optimal \code{psi} solutions composed by sigma.
#' @param H_XA A numeric vector of inverse-propensity weights, typically from \code{HX()}.
#'
#' @return A numeric vector of nu on the \[0,1\] scale.
#'
#'@export
update_nu_XA <- function(offset_nu_XA, epsilon2, sigma_psi_collection, H_XA){
  out <- expit(offset_nu_XA + H_XA*(sigma_psi_collection%*%epsilon2))
  return(out)
}

#' Update Mu via Augmented Covariate Adjustment
#'
#' Computes updated Mu predictions by adding the bias correction for all previous solutions at X 
#' \code{psi(X)} scaled by coefficients \code{epsilon1_collection},
#' then transforms back via the logistic function.
#'
#' @param A A vector of binary treatment indicators (0 or 1).
#' @param X A matrix or data frame of covariates.
#' @param mu0 A function predicting the primary outcome given treatment and covariates.
#' @param epsilon1 A numeric vector of GLM coefficients for each column in \code{psi_collection}.
#' @param theta_collection A list of the optimal \code{theta} enabling the reconstruction of optimal \code{psi} functions.
#' @param prop_score A function of A and X to compute the propensity score.
#'
#' @return A numeric vector of updated mu on the \[0,1\] scale.
#'
#'@export
update_mu <- function(A, X, mu0, epsilon1, theta_collection, prop_score){ 
  H_XA <- HX(A, X, prop_score)
  res<- lapply(theta_collection, function(theta){make_psi(theta)(X)})
  psi_collection <- do.call(cbind, res)
  out <- expit(logit(mu0(A, X))+ H_XA*(psi_collection%*%epsilon1)) 
  return(out)
}

#' Update Nu via Augmented Covariate Adjustment
#'
#' Computes updated Nu predictions by adding the bias correction for all previous solutions at X 
#' \code{sigma_beta(psi(X),beta,centered)} scaled by coefficients \code{epsilon2_collection},
#' then transforms back via the logistic function.
#'
#' @param A A vector of binary treatment indicators (0 or 1).
#' @param X A matrix or data frame of covariates.
#' @param nu0 A function predicting the adverse event given treatment and covariates.
#' @param epsilon2 A numeric vector of GLM coefficients for each column in \code{sigma_psi_collection}.
#' @param theta_collection A list of the optimal \code{theta} enabling the reconstruction of optimal \code{sigma_beta(psi)} functions.
#' @param prop_score A function of treatment and covariates to compute the propensity score.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#'
#' @return A numeric vector of updated mu on the \[0,1\] scale.
#'
#'@export
update_nu <- function(A, X, nu0, epsilon2, theta_collection, prop_score, beta=0.05, centered=FALSE){
  H_XA <- HX(A, X, prop_score)
  res<- lapply(theta_collection, function(theta){make_psi(theta)(X)})
  sigma_psi_collection <- do.call(cbind, lapply(res,function(X)sigma_beta(X, beta, centered)))
  out <- expit(logit(nu0(A,X))+ H_XA*(sigma_psi_collection%*%epsilon2))
  return(out)
}

#' Iterative Optimization Procedure
#' 
#' This function performs an iterative optimization routine to correct and minimize the objective function. 
#' It iteratively finds a solution and corrects the objective function for such optimal solution, until 
#' two consecutive solutions do not change much. 
#'
#' @param mu0 A function that returns the conditional mean of the outcome given treatment and covariates.
#' @param nu0 A function that returns the conditional mean of the adverse event.
#' @param prop_score A function that estimates the propensity score given treatment and covariates.
#' @param X A matrix or data frame of covariates.
#' @param A A binary vector (0/1) indicating treatment assignment.
#' @param Y A numeric vector or matrix of primary outcomes.
#' @param Xi A numeric vector or matrix of primary outcomes.
#' @param lambda A numeric scalar controlling the weight of the constraint function in the objective.
#' @param alpha A numeric scalar representing the constraint tolerance.
#' @param precision A numeric scalar that determines the convergence precision desired.
#' @param beta A numeric scalar controlling the sharpness of the probability function.
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta}.
#' @param filepath_results Path to the folder where intermediate and final results are saved.
#'
#' @return A list containing:
#' \item{iter}{The number of completed iterations.}
#' \item{offset_mu}{Initial logit-transformed outcome predictions.}
#' \item{offset_nu}{Initial logit-transformed auxiliary predictions.}
#' \item{psi_collection}{Matrix of covariate projections across iterations.}
#' \item{sigma_psi_collection}{Matrix of transformed projections across iterations.}
#' \item{epsilon1}{GLM coefficients from the outcome model.}
#' \item{epsilon2}{GLM coefficients from the auxiliary model.}
#' \item{theta_collection}{List of parameter vectors from each iteration of the functional weight estimation.}
#'
#' @details
#' This function saves intermediate results to files in order to recover progress or inspect iteration-level behavior. 
#' If the optimization converges or the maximum number of iterations is reached, the final parameter vector \code{theta_init} is saved.
#'
#' @examples
#' # (Requires user-defined functions: mu0, nu0, prop_score, FW, make_psi, sigma_beta, update_mu_XA, update_nu_XA)
#' # Optimization_Estimation(mu0, nu0, prop_score, df, lambda=0.1, alpha=0.1, precision=1e-4,
#' #                         beta=1, centered=TRUE, folder="path/to/folder", prefix="run1")
#'
#' @export
Optimization_Estimation <- function(mu0, nu0, prop_score, X, A, Y, Xi, lambda, alpha, precision, beta=0.05, centered=FALSE, filepath_results){
  tol <- 7.5*1e-2
  max_iter <- 1.5*1e1
  Delta_mu <- function(X){mu0(rep(1,nrow(X)),X)-mu0(rep(0,nrow(X)),X)}
  Delta_nu <- function(X){nu0(rep(1,nrow(X)),X)-nu0(rep(0,nrow(X)),X)}

  H <- HX(A,X,prop_score)
  
  theta <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
  psi<- make_psi(theta)
  psi_X <- psi(X)
  sigma_psi_X <- sigma_beta(psi_X,beta, centered)
  
  offset_mu <- qlogis(mu0(A,X))
  df_mu <- tibble::tibble(
    Y = Y)
  
  offset_nu <- qlogis(nu0(A,X))
  df_nu <- tibble::tibble(
    xi = Xi)
  
  psi_collection <- NULL
  sigma_psi_collection <- NULL
  theta_collection<-list()
  
  go_on <- TRUE
  k <- 0 
  while(go_on){
    k <- k + 1
    varname <- paste0("newcov.", k)
    df_mu[[varname]] <- H*as.vector(psi_X)
    df_nu[[varname]] <- H*as.vector(sigma_psi_X)
    psi_collection <- cbind(psi_collection, as.vector(psi_X))
    sigma_psi_collection <- cbind(sigma_psi_collection, as.vector(sigma_psi_X))
    theta_collection[[k]] <- theta
    
    mu_update_obj <- stats::glm(Y ~ -1 + ., offset=offset_mu, data = df_mu, family=binomial())
    nu_update_obj <- stats::glm(xi ~ -1 + ., offset=offset_nu, data = df_nu, family=binomial())
    
    epsilon1 <- as.matrix(as.numeric(mu_update_obj$coefficients))
    epsilon2 <- as.matrix(as.numeric(nu_update_obj$coefficients))

    Delta_mu <- function(X) { update_mu_XA(qlogis(mu0(rep(1,nrow(X)),X)), epsilon1, psi_collection, HX(rep(1,nrow(X)),X,prop_score)) - 
        update_mu_XA(qlogis(mu0(rep(0,nrow(X)),X)), epsilon1, psi_collection, HX(rep(0,nrow(X)), X, prop_score)) }
    Delta_nu <- function(X) { update_nu_XA(qlogis(nu0(rep(1,nrow(X)),X)), epsilon2, sigma_psi_collection,HX(rep(1,nrow(X)),X,prop_score)) - 
        update_nu_XA(qlogis(nu0(rep(0,nrow(X)),X)), epsilon2, sigma_psi_collection, HX(rep(0,nrow(X)), X, prop_score)) }
    
    theta <- FW(X, Delta_mu, Delta_nu, lambda, alpha, beta, centered, precision, verbose=TRUE)
    psi<- make_psi(theta)
    new_psi <- psi(X)
    sigma_psi_X <- sigma_beta(new_psi,beta, centered)
    go_on <- (k < max_iter) & (sqrt(mean((psi_X - new_psi)^2)) > tol)
    
    if(k%%10==0){
      print(mean(H*(-2*psi_X*(df_mu$Y-update_mu_XA(offset_mu, epsilon1, psi_collection, H))
                    + lambda*sigma_psi_X*(df_nu$xi - update_nu_XA(offset_nu, epsilon2, sigma_psi_collection, H)))))
      print(sqrt(mean((psi_X - new_psi)^2)))}
    
    psi_X <- new_psi
    
    out <- list(
      iter=k,
      offset_mu=offset_mu,
      offset_nu=offset_nu, 
      psi_collection=psi_collection, 
      sigma_psi_collection=sigma_psi_collection, 
      epsilon1=epsilon1, 
      epsilon2=epsilon2,
      theta_collection=theta_collection
    )
    
    step_file_prev <- file.path(paste0(filepath_results,"_step_", k - 1, ".rds"))
    step_file_current <- file.path(paste0(filepath_results,"_step_", k, ".rds"))
    
    saveRDS(out, file = step_file_current)
    
    if (file.exists(step_file_prev)) {
      file.remove(step_file_prev)
      cat("Deleted previous step file:", step_file_prev, "\n")}
  }
  return(out)
}
