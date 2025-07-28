#' Main algorithm
#' 
#' This function runs the whole algorithm. It first checks the data, then computes the nuisance parameters using SuperLearner librarie. 
#' Then proceed onto the 
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param A A binary vector of size n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param Xi A numeric vector or matrix of length n indicating adverse events (0 or 1).
#' @param Lambdas A sequence of non-negative numeric scalars controlling the penalty for violating the constraint (seq(1,8,by=1) by default).
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param precision A numeric scalar defining the desired convergence precision (0.05 by default). The number of Frank-Wolfe iterations (K) is inversely proportional to this value, calculated as 1/precision.
#' @param B A vector of non-negative scalars controlling the sharpness of the treatment probability function (c(0.05, 0.1, 0.25, 0.5) by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param root.path Path to the folder where all results are to be saved.
#' @export
Naive_approach_algorithm <- function(X, A, Y, Xi, folds, 
                                     mu0_train, mu0_test, 
                                     nu0_train, nu0_test, 
                                     prop_score_train, prop_score_test, 
                                     Lambdas=seq(1, 8, by=1), 
                                     alpha=0.1, precision=0.05,B=c(0.05, 0.1, 0.25, 0.5), 
                                     centered=FALSE, root.path){
  
  n <- nrow(X)
  # Folds for cross-validation 
  checks<- PLUCR::check_data(Y, Xi, A, X, folds) 
  ########################################
  # Step 1: Compute nuisance parameters #
  ########################################
  
  ### Training data 
  X_train <- X[folds[[2]],]
  A_train <- A[folds[[2]]]
  Y_train <- Y[folds[[2]]]
  Xi_train <- Xi[folds[[2]]]
  
  ###  Testing data 
  X_test <- X[folds[[3]],]
  A_test <- A[folds[[3]]]
  Y_test <- Y[folds[[3]]]
  Xi_test <- Xi[folds[[3]]]
  
  ##########################################
  # Step 2: Training process: grid search  #
  ##########################################
  ##### 1.1- Start by testing \code{lambda}=0 
  saved <- FALSE
  combinations <- NULL
  naive_opt_res <- NULL
  naive_opt_theta <- NULL
  stopped <- FALSE
  lambda <- 0
  min_constraint_lambda0 <- Inf
  max_policy_value <- -Inf
  beta_0 <- NULL
  
  theta_0 <- PLUCR::FW(X_train,delta_Mu=function(X){mu0_train(rep(1,nrow(X)),X)-mu0_train(rep(0,nrow(X)),X)}, 
                       delta_Nu=function(X){nu0_train(rep(1,nrow(X)),X)-nu0_train(rep(0,nrow(X)),X)}, 
                       lambda=0, beta=0.05, precision=precision) 
  
  ##### 1.2- Check whether not considering your constraint satisfies already your constraint  
  attr(theta_0, "lambda") <- 0 
  for(beta in B){
    saved <- FALSE
    res <- process_results(theta_n, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, 0, alpha,  beta, centered)
    if (!saved && res[[1]]$constraint < 0) {
      combinations <- rbind(combinations, c(beta, 0))
      naive_opt_res <- rbind(naive_opt_res,res[[1]])
      naive_opt_theta[[length(naive_opt_theta)+1]] <- theta_n
      saved <- TRUE
    }
    if(res[[2]]<0){
      stopped <- TRUE
      break
    }
  }
  attr(theta_0, "beta") <- beta_0
  ##### If your constraint was not satified with lambda=0 goes onto step 2
  ##### 2.1- Test different lambda and beta combinations and save the optimal solutions satisfying the constraint 
  ### Training 
  if(!stopped){
  for (beta in B){
    saved <- FALSE
    for (lambda in Lambdas){
      # Policy optimization
      theta_opt <- FW(X_train,delta_Mu=function(X){mu0_train(rep(1,nrow(X)),X)-mu0_train(rep(0,nrow(X)),X)}, 
                      delta_Nu=function(X){nu0_train(rep(1,nrow(X)),X)-nu0_train(rep(0,nrow(X)),X)}, 
                      lambda=lambda, beta=beta, precision=precision) 
    
      res <- process_results(theta_opt, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda, alpha,  beta, centered)
      if (!saved && res[[1]]$constraint < 0) {
        saved <- TRUE
        combinations <- rbind(combinations, c(beta, lambda))
        naive_opt_res <- rbind(naive_opt_res,res[[1]])
        naive_opt_theta[[length(naive_opt_theta)+1]] <- theta_n
        saved <- TRUE
      }
      if(res[[2]]<0){
        break
      }
    }
  }
  }
  optimal_combination <- get_opt_beta_lambda(combinations,root.path)
  files_del <- file.path(root.path,"Intermediate")
  unlink(files_del, recursive = TRUE)
  theta_final <- readRDS(file = file.path(root.path, "Theta_opt", paste0(optimal_combination$beta, "_", optimal_combination$lambda, ".rds"))) 
  
  attr(theta_final, "lambda") <- optimal_combination$lambda 
  attr(theta_final, "beta") <- optimal_combination$beta
  
  psi_values <- make_psi(theta_final)(X_test)
  optimal_treatment_rule <- sigma_beta(psi_values, beta = attr(theta_final, "beta"))
  
  # Calculate proportion over 0.5
  prop_over_0.50 <- mean(optimal_treatment_rule > 0.5)
  if(prop_0.50<0.1){
    warning(sprintf(
      paste(
        "Only %.1f%% of the test set has an optimal treatment probability above 0.5.",
        "This may indicate that your tolerance for adverse events (alpha) is too strict.",
        "Consider relaxing it if treatment is being under-assigned."), 
      100 * prop_over_0.50))
  }
  
  return(list(theta_0,theta_final))
}