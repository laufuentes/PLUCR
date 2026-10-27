expit <- stats::plogis
logit <- stats::qlogis

#' Baseline effect on Y component function.
#'
#' Computes a covariate-dependent baseline effect.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#'
#' @return A numeric vector with the transformed values based on covariates.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' b(X)
#' @export
b   <- function(X){
  out <- -tau(X) / 2 + 0.5 * (X[, 3] - 0.5)
  return(out)
} 

#' Treatment effect on Y component function.
#'
#' Computes the treatment effect component.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#'
#' @return A numeric vector with the transformed values based on covariates.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' tau(X)
#' @export
tau <- function(X){
  tau_slope = 7
  tau_vars = 1
  out <- tau_slope * (rowMeans(X[, tau_vars, drop = FALSE]) - 0.5)
  return(out)
}

#' Conditional mean of Y function.
#'
#' Computes the conditional mean of Y based on covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#' @param A A vector indicating treatment assignment (1 or 0) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates 
#' and treatment.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' A <- rep(1, 10)
#' model_Y(X, A)
#' @export
model_Y <- function(X,A){
  out <- expit(b(X) + A * tau(X))
  return(out)
}


#' Treatment effect on Xi component function.
#'
#' Computes the treatment effect component.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#'
#' @return A numeric vector with the transformed values based on covariates.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' p(X)
#' @export
p <- function(X){
  p_slope = 0.50
  out <- p_slope * X[, 2]
  return(out)
}

#'  Conditional mean of Xi function.
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#' @param A A vector indicating treatment assignment (1 or 0) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' A <- rep(1,10)
#' model_Xi(X, A)
#' @export
model_Xi <- function(X, A){
  nu_base = 0.05
  out <- nu_base + A * (1 - nu_base) * p(X)
  return(out)
}

#' Conditional Average Treatment Effect for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' delta_mu(X)
#' @export
delta_mu <- function(X){
  out <- model_Y(X, 1) - model_Y(X, 0)
  return(out)
}
attr(delta_mu, "vars")<- c(1, 3)

#' Conditional Average Treatment Effect for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(stats::runif(10*5), 10, 5)
#' delta_nu(X)
#' @export
delta_nu <- function(X){
  out <- model_Xi(X, 1) - model_Xi(X, 0)
  return(out)
}
attr(delta_nu, "vars")<- c(1, 2)

#' Synthetic data generator and functions generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param is_RCT Logical value indicating whether the scenario is an RCT (FALSE by default). 
#' @param seed Integer or NA (NA by default).
#'
#' @return A list containing two data frames (\code{df_complete} with all potential outcomes and 
#' treatment assignments and  \code{df_obs} with observed outcomes based on treatment) and the oracular 
#' functions delta_Mu and delta_Nu. 
#' @examples
#' data <- generate_data(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
generate_data <- function(n, ncov=10L, is_RCT=FALSE, seed=NA){
  ncov <- R.utils::Arguments$getIntegers(ncov, c(5, 15))
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  X <- matrix(stats::runif(n*ncov,0,1),n,ncov)
  delta_Mu <- delta_mu
  mod_Y <- model_Y
  e_slope = 2.5
  p.s <- expit(e_slope * (X[, 2] - 0.5))   
  
  if(is_RCT){
    p.s <- rep(0.5, nrow(X))
  }
  
  Treatment <- stats::rbinom(n,1,p.s)
  
  Y.1 <- stats::rbinom(n, 1, mod_Y(X,1)) 
  Y.0 <- stats::rbinom(n, 1, mod_Y(X,0))
  
  delta_Nu <- delta_nu
  mod_Xi <- model_Xi
  
  Xi.0 <- stats::rbinom(n, 1, mod_Xi(X,0))
  Xi.1<- stats::rbinom(n, 1, mod_Xi(X,1))
  
  Y_obs <- ifelse(Treatment==1, Y.1, Y.0)
  Xi_obs <- ifelse(Treatment==1, Xi.1, Xi.0)
  df_complete <- data.frame(X = X,
                            Treatment,
                            Y.1 = Y.1,
                            Y.0 = Y.0,
                            Xi.1 = Xi.1,
                            Xi.0 = Xi.0, 
                            Y = Y_obs, 
                            Xi = Xi_obs)
  df_obs<- data.frame(X = X, Treatment,
                      Y = Y_obs, Xi = Xi_obs)
  return(list(df_complete, df_obs, delta_Mu, delta_Nu))
}