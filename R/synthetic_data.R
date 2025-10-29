expit <- plogis
logit <- qlogis

#' Linear treatment effect on Y component function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_linear(X, A)
#' @export
model_Y_linear <- function(X,A){
  out <- 2*A*(1- X[,1]-X[,2])
  return(out)
}

#' Thresholded treatment effect on Y component function
#'
#' Computes a thresholded-shaped interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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
  out <- A*form
  return(out)
}

#' Mixed treatment effect on Y component function
#'
#' Computes a mix of a linear and threshold shaped interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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
  out <- -2 * A * effect
  return(out)
}

#' No treatment effect on Y component function
#'
#' Computes a null treatment effect for everyone 
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (-1 or 1).
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_null_TE(X, A)
#' @export
model_Y_null<- function(X,A){
  return(0)
}

#' Constant treatment effect on Y component function
#'
#' Computes a null treatment effect for everyone 
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (-1 or 1).
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_constant_TE(X, A)
#' @export
model_Y_constant<- function(X, A) {
  out <- A*2
  return(out)
}

#' Linear treatment effect on Xi Component Function
#'
#' Computes a linear interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_linear(X, A)
#' @export
model_Xi_linear <- function(X){
  n <- nrow(X)
  Xi.0 <- stats::rbinom(n,1,0.25)
  p1 <- expit(4*(X[,2]-1/2))
  Xi.1<- ifelse(Xi.0 == 1, 1, rbinom(n, 1, p1))
  
  return(list(Xi.0,Xi.1))
}

#' Thresholded treatment effect on Xi component function
#'
#' Computes a threshold-based interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_threshold(X, A)
#' @export
model_Xi_threshold <- function(X){
  n <- nrow(X)
  Xi.0 <- stats::rbinom(n,1,0.1)
  in_square <- (X[,3] > 0.2 & X[,3] < 0.8) & (X[,4] > 0.25 & X[,4] < 0.75)
  p1 <- ifelse(in_square, 0.85, 0.35)
  Xi.1<- ifelse(Xi.0 == 1, 1, stats::rbinom(n,1,p1))
  return(list(Xi.0,Xi.1))
}

#' Mixed treatment effect on Xi component function
#'
#' Computes a threshold-based interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_mix(X, A)
#' @export
model_Xi_mix <- function(X) {
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
  
  return(list(Xi.0,Xi.1))
}
#' Low treatment effect on Xi
#'
#' Computes a close to zero interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Xi_satisfied(X, A)
#' @export
model_Xi_satisfied <- function(X){
  n <- nrow(X)
  Xi.0 <- stats::rbinom(n,1,1e-2)
  p1 <- 4*1e-2
  Xi.1<- ifelse(Xi.0 == 1, 1, rbinom(n, 1, p1))
  return(list(Xi.0,Xi.1))
}

#' Linear-shaped Conditional Average Treatment Effect estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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

#' Thresholded-shaped Conditional Average Treatment Effect estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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

#' Mixed-shape Conditional Average Treatment Effect estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control, using \code{h_Y}.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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

#' Null Conditional Average Treatment Effect estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu_null(X)
#' @export
delta_mu_null <- function(X){
  n <- nrow(X)
  out <- 0.55*(expit(model_Y_null(X,rep(1,n)))-expit(model_Y_null(X,rep(-1,n))))
  return(out)
}
attr(delta_mu_null, "vars")<- c(1, 2)

#' Constant Conditional Average Treatment Effect estimator for Y
#'
#' Computes the difference in expected Y outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#'
#' @return A numeric vector that represents the contrast between primary outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_mu_null(X)
#' @export
delta_mu_constant <- function(X){
  n <- nrow(X)
  out <- 0.55*(expit(model_Y_constant(X,rep(1,n)))-expit(model_Y_constant(X,rep(-1,n))))
  return(out)
}
attr(delta_mu_constant, "vars")<- c(1, 2)


#' Linear-shaped Conditional Average Treatment Effect estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu(X)
#' @export
delta_nu_linear <- function(X){
  p0 <- 0.25
  p1 <- 0.25 + 0.75*expit(4*(X[,2]-1/2)) 
  out <- p1-p0
  return(out)
}
attr(delta_nu_linear, "vars")<- c(1, 2)

#' Thresholded Conditional Average Treatment Effect estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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
  out <- p1-p0
  return(out)
}
attr(delta_nu_threshold, "vars")<- c(3, 4)

#' Mixed-shaped Conditional Average Treatment Effect estimator for Xi
#'
#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
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

#' Computes the difference in expected outcomes under treatment and control.
#'
#' @param X A matrix of covariates of size n x d (input data in [0,1]).
#'
#' @return A numeric vector that represents the contrast between adverse event outcomes for given \code{X}.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' delta_nu_mix(X)
#' @export
delta_nu_satisfied <- function(X){
  p0 <- 1e-2
  effect <- 4*1e-2
  p1 <- p0 + (1-p0)*effect
  out <- p1-p0
  return(out)
}
attr(delta_nu_satisfied, "vars")<- c(1, 2)

#' Synthetic data generator and functions generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param scenario_mu String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix", "Null", "Constant").
#' @param scenario_nu String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix", "Satisfied").
#' @param is_RCT Logical value indicating whether the scenario is an RCT (FALSE by default). 
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
generate_data <- function(n, ncov=10L, scenario_mu=c("Linear", "Threshold", "Mix", "Null", "Constant"), 
                          scenario_nu=c("Linear", "Threshold", "Mix", "Satisfied"), is_RCT=FALSE, 
                          seed=NA){
  # ncov <- 5L  
  ncov <- R.utils::Arguments$getIntegers(ncov, c(5, 15))
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
  }else if(scenario_mu=="Mix"){
    delta_Mu <- delta_mu_mix
    mod_Y <- model_Y_mix 
    p.s <- expit(0.5*X[,3]+0.5*ifelse(X[,1]>0.3,1,0))
  } else if(scenario_mu=="Null"){
    delta_Mu <- delta_mu_null
    mod_Y <- model_Y_null 
    p.s <- expit(4*(X[,5]-1/2))
  } else{
    delta_Mu <- delta_mu_constant
    mod_Y <- model_Y_constant
    p.s <- expit(4*(X[,5]-1/2))
  }
  
  if(is_RCT){
    p.s <- rep(0.5, nrow(X))
  }
  outcome_X <- 3*X[,3] - X[,4]
  epsilon_Y <- rnorm(n,0,1)  
  Treatment <- stats::rbinom(n,1,p.s)
  
  if(scenario_mu %in% c("Linear", "Threshold", "Mix")){
    Y.1 <- 0.05 * expit(epsilon_Y) + 
      0.95 * expit(mod_Y(X,rep(1,n)))
    Y.0 <- 0.05 * expit(epsilon_Y) + 
      0.95 * expit(mod_Y(X,rep(-1,n)))
  }else{
    Y.1 <- 0.05 * expit(epsilon_Y) + 0.55*expit(mod_Y(X,rep(1,n))) +
      0.35*expit(outcome_X) 
    Y.0 <- 0.05 * expit(epsilon_Y) + 0.55*expit(mod_Y(X,rep(-1,n))) +
      0.35*expit(outcome_X)
  }
  if(scenario_nu=="Linear"){
    delta_Nu <- delta_nu_linear
    mod_Xi <- model_Xi_linear
  }else if(scenario_nu=="Threshold"){
    delta_Nu <- delta_nu_threshold
    mod_Xi <- model_Xi_threshold
  }else if(scenario_nu=="Mix"){
    delta_Nu <- delta_nu_mix
    mod_Xi <- model_Xi_mix 
  }else{
    delta_Nu <- delta_nu_satisfied
    mod_Xi <- model_Xi_satisfied  
  }
  res_xi <- mod_Xi(X)
  Xi.0 <- res_xi[[1]]
  Xi.1<- res_xi[[2]]
  
  df_complete <- data.frame(X=X,
                            Treatment,
                            Y.1=Y.1,
                            Y.0=Y.0,
                            Xi.1=Xi.1,
                            Xi.0=Xi.0, 
                            Y=ifelse(Treatment==1,Y.1,Y.0), 
                            Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  df_obs<- data.frame(X=X,
                      Treatment,
                      Y=ifelse(Treatment==1,Y.1,Y.0),
                      Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  return(list(df_complete, df_obs, delta_Mu, delta_Nu))
}

#' Realistic treatment effect on Y component function
#'
#' Computes a realistic interaction term between covariates and treatment.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A vector indicating treatment assignment (+1 or -1) for each observation.
#'
#' @return A numeric vector with the transformed values based on covariates and treatment.
#' @examples
#' X <- matrix(runif(10*2), 10, 2)
#' A <- rep(1, 10)
#' model_Y_linear(X, A)
#' @export
model_Y_realistic <- function(X,A){
  out <- A * (-4 + 0.1 * X[,1])
  return(out)
}

#' Realistic treatment effect on Xi Component Function
#'
#' Computes a realistic interaction term between covariates and treatment.
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
model_Xi_realistic <- function(X){
  n <- nrow(X)
  p0 <- 0.01 
  xi.0 <- rbinom(n,1,p0)
  p1 <- ifelse(X[,2]==0, p0 , 0.35)
  p1 <- ifelse((X[,3] == 1),1,p1)
  xi.1 <- rbinom(n,1,p1)
  xi.1 <- ifelse(xi.0==1, 1, xi.1)
  return(list(xi.0, xi.1))
}

#' Realistic Conditional Average Treatment Effect estimator for Y
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
delta_mu_realistic <- function(X){
  n <- nrow(X)
  temp <- (model_Y_realistic(X,rep(1,n)))-(model_Y_realistic(X,rep(-1,n)))
  out <- temp / (attr(X,"max_Y")- attr(X,"min_Y"))
  return(out)
}
attr(delta_mu_realistic, "vars")<- c(1, 4)

#' Realistic Conditional Average Treatment Effect estimator for Xi
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
delta_nu_realistic <- function(X){
  p0 <- 0.01
  p1 <- ifelse(X[,3]==1,1, ifelse(X[,2]==1,0.35,0))
  out <- p1-p0
  return(out)
}
attr(delta_nu_realistic, "vars")<- c(2, 3)


#' Realistic synthetic data generator and functions generator
#'
#' Generates a realistic dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param ncov Number of baseline covariates (at least 2L and 10L by default).
#' @param scenario_mu String indicating the type of scenario for delta_Mu ("Linear", "Threshold", "Mix", "Null", "Constant").
#' @param scenario_nu String indicating the type of scenario for delta_Nu ("Linear", "Threshold", "Mix", "Satisfied").
#' @param is_RCT Logical value indicating whether the scenario is an RCT (FALSE by default). 
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
generate_realistic_data <- function(n, ncov=5L, scenario_mu="Realistic", scenario_nu="Realistic", 
                                    is_RCT=FALSE, seed=NA){
  ncov <- R.utils::Arguments$getIntegers(ncov, c(2, 15))
  scenario_mu <- match.arg(scenario_mu)
  scenario_nu <- match.arg(scenario_nu)
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  sex <- as.integer(rbinom(n,1, 0.5)) # 1: woman, 0: man
  age <- round(runif(n,16,65))
  is_pregnancy_window <- ifelse(age >= 18 & age <= 45 & sex == 1, 1, 0)
  is_pregnant <- ifelse(is_pregnancy_window==0, 0, rbinom(n, 1, 0.3))
  
  X <- matrix(cbind(age, 
                    sex, 
                    is_pregnant,
                    runif(n, min = 0, max=10), 
                    runif(n, min = -5, max=5)), 
              n,ncov)
  colnames(X) <- c("1", "2", "3", "4", "5")
  delta_Mu <- delta_mu_realistic
  delta_Nu <- delta_nu_realistic
  mod_Y <- model_Y_realistic
  mod_Xi <- model_Xi_realistic
  p.s <- expit(-0.5*sex + 0.2 * (X[,5]) + 0.6 * (X[,4]-5.5))
  
  if(is_RCT){
    p.s <- rep(0.5, nrow(X))
  }
  outcome_X <- (0.4*X[,4] - 0.2*X[,5])
  epsilon_Y <- rnorm(n,0,1)  
  Treatment <- stats::rbinom(n,1,p.s)
  
  Y.1 <- 0.5*epsilon_Y + mod_Y(X,rep(1,n)) + outcome_X
  Y.0 <- 0.5*epsilon_Y + mod_Y(X,rep(-1,n)) + outcome_X
  delta_Mu <- delta_mu_realistic
  delta_Nu <- delta_nu_realistic
  
  res_xi <-model_Xi_realistic(X)
  Xi.0 <- res_xi[[1]]
  Xi.1 <- res_xi[[2]]
  
  df_complete <- data.frame(X=X,Treatment, Y.1=Y.1, Y.0=Y.0,
                            min_Y=min(c(Y.1,Y.0)), 
                            max_Y=max(c(Y.1,Y.0)),
                            Xi.1=Xi.1, Xi.0=Xi.0, 
                            Y=ifelse(Treatment==1,Y.1,Y.0), 
                            Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  df_obs<- data.frame(X=X,Treatment,
                      Y=ifelse(Treatment==1,Y.1,Y.0),
                      Xi=ifelse(Treatment==1,Xi.1,Xi.0))
  
  return(list(df_complete, df_obs, delta_Mu, delta_Nu))
}