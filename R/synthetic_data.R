expit <- plogis
logit <- qlogis

#' model_Y_linear: Treatment Effect on Y Component Function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' h_Y(X, A)
#' @export
model_Y_linear <- function(X,A){
  return(2*A*(1- X[,1]-X[,2]))
}

#' model_Y_threshold: Thresholded Treatment Effect on Y Component Function
#'
#' Computes a thresholded-shaped interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector or matrix of length n indicating treatment assignment (-1 or 1).
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_threshold(X, A)
#' @export
model_Y_threshold <- function(X, A) {
  form <- ifelse(X[,1]>0.4,0.6,-0.4) + ifelse(X[,2]>0.6, 0.5,-0.5)
  return(A*form)
}

#' model_Y_mix: Mixed Treatment Effect on Y Component Function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector or matrix of length n indicating treatment assignment (-1 or 1).
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_mix(X, A)
#' @export
model_Y_mix <- function(X, A) {
  threshold <- 0.3
  linear_pred <- -(1- X[,1]-X[,2])
  effect <- ifelse(X[,1] < threshold & X[,2] < threshold,-3, linear_pred)
  return(-2 * A * effect)
}

#' model_Xi_linear: Treatment Effect on Xi Component Function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_linear(X, A)
#' @export
model_Xi_linear <- function(X,A){
  n <- nrow(X)
  Xi.0 <- stats::rbinom(n,1,0.25)
  p1 <- expit(4*(X[,2]-1/2))
  Xi.1<- ifelse(Xi.0 == 1, 1, rbinom(n, 1, p1))
  return(ifelse(A==1, Xi.1, Xi.0))
}

#' model_Xi_threshold: Treatment Effect on Xi Component Function
#'
#' Computes a threshold-based interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_linear(X, A)
#' @export
model_Xi_threshold <- function(X,A){
  n <- nrow(X)
  Xi.0 <- stats::rbinom(n,1,0.1)
  in_square <- (X[,3] > 0.2 & X[,3] < 0.8) & (X[,4] > 0.25 & X[,4] < 0.75)
  p1 <- ifelse(in_square, 0.85, 0.35)
  Xi.1<- ifelse(Xi.0 == 1, 1, stats::rbinom(n,1,p1))
  return(ifelse(A==1, Xi.1, Xi.0))
}

#' model_Xi_threshold: Treatment Effect on Xi Component Function
#'
#' Computes a threshold-based interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_linear(X, A)
#' @export
model_Xi_mix <- function(X, A) {
  n <- nrow(X)
  threshold <- 0.3
  
  # Distance from center
  d <- sqrt((X[,1] - 0.5)^2 + (X[,2] - 0.5)^2)
  
  # Linear predictor outside the circle
  linear_pred <- 4*(X[,2]-1/2)
  
  # Strong effect inside the circular region
  effect <- ifelse(d < threshold, -3, linear_pred)
  
  # Convert to probability with expit (to stay in [0,1])
  prob <- expit(effect)
  
  # Simulate binary potential outcomes
  Xi.0 <- rbinom(n, 1, 0.1)
  Xi.1 <- rbinom(n, 1, prob)
  
  return(ifelse(A == 1, Xi.1, Xi.0))
}

#' delta_mu_linear: Conditional Average Treatment Effect Estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu(X)
#' @export
delta_mu_linear <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(model_Y_linear(X,rep(1,n)))-expit(model_Y_linear(X,rep(-1,n))))
  return(out)
}
attr(delta_mu_linear, "vars")<- c(1, 2)

#' delta_mu_threhsold: Thresholded-shaped Conditional Average Treatment Effect Estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu(X)
#' @export
delta_mu_threshold <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(model_Y_threshold(X,rep(1,n)))-expit(model_Y_threshold(X,rep(-1,n))))
  return(out)
}
attr(delta_mu_threshold, "vars")<- c(1, 2)

#' delta_mu_mix: Mixed conditional Average Treatment Effect Estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu(X)
#' @export
delta_mu_mix <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(model_Y_mix(X,rep(1,n)))-expit(model_Y_mix(X,rep(-1,n))))
  return(out)
}
attr(delta_mu_mix, "vars")<- c(1, 2)

#' delta_nu_linear: Conditional Average Treatment Effect Estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu(X)
#' @export
delta_nu_linear <- function(X){
  p0 <- 0.25
  p1 <- 0.25 + 0.75*expit(4*(X[,2]-1/2))
  return(p1-p0)
}
attr(delta_nu_linear, "vars")<- c(1, 2)

#' delta_nu_threshold: Conditional Average Treatment Effect Estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu(X)
#' @export
delta_nu_threshold <- function(X){
  p0 <- 0.1
  in_square <- (X[,3] > 0.2 & X[,3] < 0.8) & (X[,4] > 0.25 & X[,4] < 0.75)
  p1 <- 0.1+  0.9*ifelse(in_square, 0.85, 0.35)
  return(p1-p0)
}
attr(delta_nu_threshold, "vars")<- c(3, 4)

#' delta_nu_mix: Mixed Conditional Average Treatment Effect Estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu_mix(X)
#' @export
delta_nu_mix <- function(X){
  p0 <- 0.1
  threshold <- 0.3
  d <- sqrt((X[,1] - 0.5)^2 + (X[,2] - 0.5)^2)
  linear_pred <- 4*(X[,2]-1/2)
  effect <- ifelse(d < threshold, -3, linear_pred)
  p1 <- p0 + (1-p0)*expit(effect)
  return(p1-p0)
}
attr(delta_nu_mix, "vars")<- c(1, 2)


#' generate_data: Synthetic Data Generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param scenario_mu String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix").
#' @param scenario_nu String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix").
#' @param seed Integer or NA (NA by default).
#'
#' @return A list containing two data frames (\code{df_complete} with all potential outcomes and 
#' treatment assignments and  \code{df_obs} with observed outcomes based on treatment) and the oracular 
#' functions delta_Mu and delta_Nu. 
#' @examples
#' data <- data_gen(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
generate_data <- function(n, ncov=10L, scenario_mu=c("Linear", "Threshold", "Mix"), scenario_nu=c("Linear", "Threshold", "Mix"), seed=NA){
  ncov <- R.utils::Arguments$getIntegers(ncov, c(2, 15))
  scenario_mu <- match.arg(scenario_mu)
  scenario_nu <- match.arg(scenario_nu)
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  X <- matrix(stats::runif(n*ncov,0,1),n,ncov)
  if(scenario_mu=="Linear"){
    delta_Mu <- delta_mu_linear
    mod_Y <- model_Y_linear
    p.s <- expit(4*(X[,2]-1/2))
  }else if(scenario_mu=="Threshold"){
    delta_Mu <- delta_mu_threshold
    mod_Y <- model_Y_threshold
    p.s <- expit(4*(X[,2]-1/2))
  }else{
    delta_Mu <- delta_mu_mix
    mod_Y <- model_Y_mix 
    p.s <- expit(0.5*X[,3]+0.5*ifelse(X[,1]>0.3,1,0))
  }
  Treatment <- stats::rbinom(n,1,p.s)
  epsilon_Y <- rnorm(n,0,1)
  
  Y.1 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(mod_Y(X,rep(1,n)))
  Y.0 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(mod_Y(X,rep(-1,n)))
  
  if(scenario_nu=="Linear"){
    delta_Nu <- delta_nu_linear
    mod_Xi <- model_Xi_linear
  }else if(scenario_nu=="Threshold"){
    delta_Nu <- delta_nu_threshold
    mod_Xi <- model_Xi_threshold
  }else{
    delta_Nu <- delta_nu_mix
    mod_Xi <- model_Xi_mix 
  }
  Xi.0 <- mod_Xi(X, rep(-1,n))
  Xi.1<- mod_Xi(X, rep(1,n))
  
  df_complete <- data.frame(X=X,
                            Treatment,
                            Y.1=Y.1,
                            Y.0=Y.0,
                            Xi.1=Xi.1,
                            Xi.0=Xi.0)
  df_obs<- data.frame(X=X,
                      Treatment,
                      Y=ifelse(Treatment==1,Y.1,Y.0),
                      Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs, delta_Mu, delta_Nu))
}