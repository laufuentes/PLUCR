expit <- plogis
logit <- qlogis
#' h_Y: Treatment Effect on Y Component Function
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
h_Y<- function(X,A){
  return(2*(1- X[,1]-X[,2])*A)
}

#' h_Y_complicated: Complicated Treatment Effect on Y Component Function
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
#' h_Y_complicated(X, A)
#' @export
h_Y_complicated <- function(X, A) {
  a <- -5
  c <- 0.5
  d <- 0.6
  parab_val <- X[,2] - (a * (X[,1] - c)^2 + d)
  return(-2 * A * parab_val)
}

#' h_Y_tree: Tree Treatment Effect on Y Component Function
#'
#' Computes a tree-shaped interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector or matrix of length n indicating treatment assignment (-1 or 1).
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' h_Y_tree(X, A)
#' @export
h_Y_tree <- function(X, A) {
   form <- ifelse(X[,1]>0.4,0.6,-0.4) + ifelse(X[,2]>0.6, 0.5,-0.5)
  return(A*form)
}

#' Synthetic Data Generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param seed Integer or NA (NA by default).
#'
#' @return A list containing two data frames: \code{df_complete} with all potential outcomes and 
#' treatment assignments, and \code{df_obs} with observed outcomes based on treatment.
#' @examples
#' data <- data_gen(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
data_gen <- function(n,seed=NA){
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  X <- matrix(stats::runif(n*10,0,1),n,10)
  p.A <- expit(4*(X[,2]-1/2))
  Treatment <- stats::rbinom(n,1,p.A)
  epsilon_Y <- stats::rnorm(n,0,1)
  
  Y.1 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y(X,rep(1,n)))
  Y.0 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y(X,rep(-1,n)))
  # logit_Xi0 <- -2 + 1.5 * X[, 3] - 1.2 * X[, 4]
  # p.Xi0 <- expit(logit_Xi0)
  # logit_r <- 1.5 * X[, 2] - 0.5
  # r <- expit(logit_r)
  # p.Xi1 <- p.Xi0 + (1 - p.Xi0) * r  # guarantees p.Xi1 ≥ p.Xi0
  # U <- runif(n)
  # Xi.0 <- as.integer(U < p.Xi0)
  # Xi.1 <- as.integer(U < p.Xi1)
  Xi.0 <- stats::rbinom(n,1,0.25)
  p1 <- expit(4*(X[,2]-1/2))
  Xi.1<- ifelse(Xi.0 == 1, 1, rbinom(n, 1, p1))
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,Xi.1=Xi.1,Xi.0=Xi.0)
  df_obs<- data.frame(X=X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs))
}

#' Complicated synthetic Data Generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param seed Integer or NA (NA by default).
#'
#' @return A list containing two data frames: \code{df_complete} with all potential outcomes and 
#' treatment assignments, and \code{df_obs} with observed outcomes based on treatment.
#' @examples
#' data <- data_gen_complicated(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
data_gen_complicated <- function(n,seed=NA){
  if(!is.na(seed)){
    set.seed(seed)
  }
  X <- matrix(stats::runif(n*10,0,1),n,10)
  p.A <- expit(0.5*X[,3]+0.5*ifelse(X[,1]>0.3,1,0))
  Treatment <- stats::rbinom(n,1,p.A)
  epsilon_Y <- stats::rnorm(n,0,1)
  
  Y.1 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y_complicated(X,rep(1,n)))
  Y.0 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y_complicated(X,rep(-1,n)))
  
  Xi.0 <- rbinom(n,1,0.1)
  d <- sqrt((X[,1]- 0.5)^2 + (X[,2] - 0.5)^2)
  # Treatment effect is large inside the circle (radius ~ 0.3), low outside
  p1 <- expit(20 * (d - 0.3))
  Xi.1<- ifelse(Xi.0==1,1,stats::rbinom(n,1,p1))
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,Xi.1=Xi.1,Xi.0=Xi.0)
  df_obs<- data.frame(X=X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs))
}

#' Synthetic Data Generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param seed Integer or NA (NA by default).
#'
#' @return A list containing two data frames: \code{df_complete} with all potential outcomes and 
#' treatment assignments, and \code{df_obs} with observed outcomes based on treatment.
#' @examples
#' data <- data_gen(100)
#' head(data[[1]])  # complete data
#' head(data[[2]])  # observed data
#' @export
data_gen_tree <- function(n,seed=NA){
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  X <- matrix(stats::runif(n*10,0,1),n,10)
  p.A <- expit(4*(X[,2]-1/2))
  Treatment <- stats::rbinom(n,1,p.A)
  epsilon_Y <- stats::rnorm(n,0,1)
  
  Y.1 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y_tree(X,rep(1,n)))
  Y.0 <- 0.05 * expit(epsilon_Y) + 
    0.95 * expit(h_Y_tree(X,rep(-1,n)))
  
  Xi.0 <- stats::rbinom(n,1,0.1)
  in_square <- (X[,3] > 0.2 & X[,3] < 0.8) & (X[,4] > 0.25 & X[,4] < 0.75)
  # Treatment effect is large inside the circle (radius ~ 0.3), low outside
  p1 <- ifelse(in_square, 0.85, 0.35)
  Xi.1<- ifelse(Xi.0 == 1, 1, stats::rbinom(n,1,p1))
  df_complete <- data.frame(X=X,Treatment,y1=Y.1,y0=Y.0,Xi.1=Xi.1,Xi.0=Xi.0)
  df_obs<- data.frame(X=X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs))
}

#' Complicated conditional Average Treatment Effect Estimator for Y
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
delta_mu_complicated <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(h_Y_complicated(X,rep(1,n)))-expit(h_Y_complicated(X,rep(-1,n))))
  return(out)
}

#' Conditional Average Treatment Effect Estimator for Y
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
delta_mu <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(h_Y(X,rep(1,n)))-expit(h_Y(X,rep(-1,n))))
  return(out)
}

#' Tree-shaped Conditional Average Treatment Effect Estimator for Y
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
delta_mu_tree <- function(X){
  n <- nrow(X)
  out <- 0.95*(expit(h_Y_tree(X,rep(1,n)))-expit(h_Y_tree(X,rep(-1,n))))
  return(out)
}

#' Conditional Average Treatment Effect Estimator for Xi
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
delta_nu <- function(X){
  p0 <- 0.25
  p1 <- 0.25 + 0.75*expit(4*(X[,2]-1/2))
  return(p1-p0)
}

#' Complicated Conditional Average Treatment Effect Estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu_complicated(X)
#' @export
delta_nu_complicated <- function(X){
  p0 <- 0.05
  d <- sqrt((X[,1]- 0.5)^2 + (X[,2] - 0.5)^2)
  # Treatment effect is large inside the circle (radius ~ 0.3), low outside
  p1 <- 0.1+ 0.9*expit(20 * (d - 0.3))
  return(p1-p0)
}

#' Conditional Average Treatment Effect Estimator for Xi
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
delta_nu_tree <- function(X){
  p0 <- 0.1
  in_square <- (X[,3] > 0.2 & X[,3] < 0.8) & (X[,4] > 0.25 & X[,4] < 0.75)
  p1 <- 0.1+ 0.9*ifelse(in_square, 0.85, 0.35)
  return(p1-p0)
}
