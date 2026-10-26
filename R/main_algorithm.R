utils::globalVariables(c(
  "SL.library", "theta_oracular", "X_test", "df_complete",
  "n", "metric", "value", "method", "type",
  "delta_mu", "delta_nu", "x", "y", "Values", "treat_proba"
))
#' Main algorithm
#'
#' Executes the full estimation pipeline for learning an optimal treatment policy under constraints.
#' This function performs the complete procedure implemented in the PLUCR framework. 
#' It begins by validating and preprocessing the input data, followed by estimating nuisance parameters using the `SuperLearner` ensemble library. 
#' Subsequently, it estimates the optimal treatment policy via the Frank-Wolfe optimization algorithm. 
#' The procedure includes an inner grid search over candidate values of `lambda` and `beta` to identify the policy 
#' that maximizes the expected primary outcome (policy value) while satisfying a constraint on the expected rate of adverse events.
#' 
#'
#' @param X A matrix of covariates of size n x d (input data in `[0,1]`).
#' @param A A binary vector of size n indicating treatment assignment (0 or 1).
#' @param Y A numeric vector or matrix of length n representing primary outcomes 
#' (in `[0,1]`).
#' @param Xi A numeric vector or matrix of length n indicating adverse events 
#' (0 or 1).
#' @param Lambdas A sequence of non-negative numeric scalars controlling the 
#' penalty for violating the constraint (seq(1,8,by=1) by default).
#' @param alpha A numeric scalar representing the constraint tolerance (in 
#' `[0,1/2]`, 0.1 by default).
#' @param precision A numeric scalar defining the desired convergence precision 
#' (0.05 by default). The number of Frank-Wolfe iterations (K) is inversely proportional to this value, calculated as 1/precision.
#' @param B A vector of non-negative scalars controlling the sharpness of the 
#' treatment probability function (c(0, 0.05, 0.1, 0.25, 0.5) by default).
#' @param centered A logical value indicating whether to apply centering in 
#' \code{sigma_beta} (FALSE by default).
#' @param Jfold Number of folds for the main algorithm, needs to be set to 3L.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param conf_level A numeric scalar representing the confidence level 
#' (0.05 by default).
#' @param SL.library1 Vector of libraries for training SuperLearner 
#' (c("SL.mean","SL.glm","SL.ranger","SL.grf") by default). This library is used 
#' to fit the nuisances for training PLUC.
#' @param SL.library2 Vector of libraries for training SuperLearner 
#' (c("SL.mean","SL.glm","SL.ranger","SL.grf") by default). This library is used 
#' to evaluate policies.
#' @param tol A numeric scalar used as an early stopping criterion based on the 
#' RMSE between consecutive solutions (0.025 by default).
#' @param max_iter A numeric scalar specifying the maximum number of iterations 
#' (5 by default).
#' @param root.path Path to the folder where all results are to be saved.
#' 
#' @return A list of matrices (`theta_0` and `theta_final`), or `theta_0` alone. 
#' These matrices are used to construct the optimal treatment rule in two steps. 
#' First, build `psi` using the `make_psi` function and evaluate it at `X` (i.e., `psi(X)`). 
#' Then, obtain the optimal treatment rule by applying \code{sigma_beta} to the selected `beta` attribute.
#' @export
main_algorithm <- function(X, A, Y, Xi, 
                           Lambdas=seq(1, 8, by=1), alpha=0.1, precision=0.05,
                           B=c(0, 0.05, 0.1, 0.25, 0.5), centered=FALSE,
                           Jfold=3L, V=2L, conf_level = 0.05,
                           SL.library1=c("SL.mean","SL.glm","SL.ranger","SL.grf"), 
                           SL.library2=c("SL.mean","SL.glm","SL.ranger","SL.grf"), 
                           tol=0.025, max_iter=5, root.path){
  # Check whether the root.path exists and contains proper folder to save data
  if (!dir.exists(root.path)) {
    warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
    dir.create(root.path, recursive = TRUE)
  }
  # Subdirectories to check
  subdirs <- c("Mu.hat", "Nu.hat", "PS.hat", "Folds", "Images", "Evaluation", "Theta_opt")
  
  for (subdir in subdirs) {
    subdir_path <- file.path(root.path, subdir)
    if (!dir.exists(subdir_path)) {
      warning(sprintf("Subdirectory '%s' is missing. Creating it...", subdir_path))
      dir.create(subdir_path, recursive = TRUE)
    }
  }
  message("Directory check complete.")
  
  n <- nrow(X)
  # Folds for cross-validation 
  folds <- SuperLearner::CVFolds(n, 
                                 id = NULL,
                                 Y = Y,
                                 cvControl = SuperLearner::SuperLearner.CV.control(V = Jfold, 
                                                                                   shuffle = TRUE))
  #Save folds
  saveRDS(folds, file=file.path(root.path,"Folds","folds.rds")) 
  # Check data
  checks<- PLUCR::check_data(Y, Xi, A, X, folds) 
  ########################################
  # Step 1: Compute nuisance parameters #
  ########################################
  # Primary outcome model 
  mu.hat.nj <- PLUCR::estimate_mu(Y=Y, A=A, X=X, folds=folds,
                                  SL.library1 = SL.library1, 
                                  SL.library2 = SL.library2, V=2L)
  #Save primary outcome model
  saveRDS(mu.hat.nj, file=file.path(root.path,"Mu.hat","mu.hat.nj.rds")) 
  # Adverse event model
  nu.hat.nj <- PLUCR::estimate_nu(Xi=Xi, A=A, X=X, folds=folds, 
                                  SL.library1 = SL.library1,
                                  SL.library2 = SL.library2, V=2L)
  #Save adverse events outcome model
  saveRDS(nu.hat.nj, file=file.path(root.path,"Nu.hat","nu.hat.nj.rds")) 
  # Propensity score
  ps.hat.nj <- PLUCR::estimate_ps(A=A, X=X, folds=folds, 
                                  SL.library = SL.library2, V=2L)
  #Save propensity score outcome model
  saveRDS(ps.hat.nj, file=file.path(root.path,"PS.hat","ps.hat.nj.rds")) 
  
  ### Training data & functions
  # data
  X_train <- X[folds[[2]],]
  A_train <- A[folds[[2]]]
  Y_train <- Y[folds[[2]]]
  Xi_train <- Xi[folds[[2]]]
  # functions
  mu0_train <- function(a,x){mu.hat.nj(a=a,x=x,ff=1)}
  nu0_train <- function(a,x){nu.hat.nj(a=a,x=x,ff=1)}
  prop_score_train <- function(a,x){ps.hat.nj(a=a,x=x,ff=1)}
  
  ###  Testing data & functions
  # data
  X_test <- X[folds[[3]],]
  A_test <- A[folds[[3]]]
  Y_test <- Y[folds[[3]]]
  Xi_test <- Xi[folds[[3]]]
  # functions
  mu0_test <- function(a,x){mu.hat.nj(a=a,x=x,ff=3)}
  nu0_test <- function(a,x){nu.hat.nj(a=a,x=x,ff=3)}
  prop_score_test <- function(a,x){ps.hat.nj(a=a,x=x,ff=3)}
  
  ##########################################
  # Step 2: Training process: grid search  #
  ##########################################
  first_lambda <- Lambdas[1]
  Lambdas <- sort(Lambdas, decreasing = TRUE)
  cardLambda <- length(Lambdas) + 1 # add the test of lambda = 0 
  cardB <- length(B)
  saved <- FALSE
  combinations <- NULL
  ##### 1.1- Test different \code{lambda} and \code{beta} combinations and save 
  ##### the last solution satisfying the constraint
  for (beta in B){
    last_sucessfull_lambda <- saved_result <- theta_opt <- NULL
    for(lambda in Lambdas){
      if(beta * lambda >= 8){break}
      out <- PLUCR::Optimization_Estimation(mu0 = mu0_train, nu0 = nu0_train, 
                                            prop_score = prop_score_train, 
                                            X = X_train, A = A_train, 
                                            Y = Y_train, Xi = Xi_train, 
                                            lambda = lambda, alpha = alpha, 
                                            precision = precision, beta = beta, 
                                            centered = centered, tol = tol)
      theta <- out$theta_collection[[length(out$theta_collection)]]
      res <- process_results(theta = theta_opt, X = X_test, A = A_test, 
                             Y = Y_test, Xi = Xi_test, mu0 = mu0_test, 
                             nu0 = nu0_test, prop_score = prop_score_test, 
                             lambda = lambda, alpha = alpha, beta = beta, 
                             centered = centered, conf_level = conf_level, 
                             card_lambda = cardLambda, card_B = cardB)
      if(res$upper_bound_constraint <= 0){
        # Save solutions satisfying the constraint
        saved_res <- res 
        theta_opt <- theta
        attr(theta_opt, "lambda") <- lambda
        attr(theta_opt, "beta") <- beta
        last_success_lambda <- lambda
      }else{break}
    }
    # Keep only the last solution satisfying the constraint
    if (!is.null(last_success_lambda)) {
      saveRDS(saved_result, 
              file = file.path(root.path, "Evaluation", 
                               paste0(beta, "_", last_success_lambda, ".rds")))
      saveRDS(theta_opt, 
              file = file.path(root.path, "Theta_opt", 
                               paste0(beta, "_", last_success_lambda, ".rds")))
      combinations <- rbind(combinations, c(beta, last_success_lambda))
      saved <- TRUE
    }
  }
  ##### 1.2- Test \code{lambda}=0 if all previous solutions satisfied the constraint
  if(last_success_lambda == first_lambda){
    out0 <- PLUCR::Optimization_Estimation(mu0=mu0_train, nu0=nu0_train, 
                                           prop_score=prop_score_train, 
                                           X=X_train, A=A_train, Y=Y_train, 
                                           Xi=Xi_train, lambda=0, alpha=alpha, 
                                           precision=precision, beta=0, 
                                           centered=centered, tol=tol, 
                                           max_iter=max_iter)
    theta_0 <- out0$theta_collection[[length(out0$theta_collection)]]
    attr(theta_0, "lambda") <- 0
    for (beta in B){
      res_0 <- process_results(theta= theta_0, X = X_test, A = A_test, 
                               Y = Y_test, Xi = Xi_test, mu0 = mu0_test, 
                               nu0 = nu0_test, prop_score = prop_score_test, 
                               lambda=0, alpha = alpha, beta = beta, 
                               centered = centered, conf_level = conf_level,
                               card_lambda = cardLambda, card_B = cardB)
      # Loop to check constraint satisfaction
      if (res_0$upper_bound_constraint <= 0) {
        saveRDS(res_0, file = file.path(root.path, "Evaluation", 
                                        paste0(beta, "_", 0, ".rds")))
        attr(theta_0, "beta") <- beta
        saveRDS(theta_0, file = file.path(root.path, "Theta_opt", 
                                          paste0(beta, "_", 0, ".rds")))
        combinations <- rbind(combinations, c(beta, 0))
        saved <- TRUE
      }
    }
  }
  if (!saved) {warning("Algorithm stopped: None of the parameters 
                        satisfied the constraint")}

  if(!is.null(combinations)){
    # Select the optimal combination (beta, lambda)
    optimal_combination <- get_opt_beta_lambda(combinations,root.path)
    beta <- optimal_combination$beta
    lambda <- optimal_combination$lambda
    theta_keep <- paste0(beta, "_", lambda, ".rds")
    
    # Delete unwanted files in Theta_opt
    theta_files <- list.files(file.path(root.path, "Theta_opt"), pattern="^\\d+")
    theta_to_delete <- theta_files[basename(theta_files) != theta_keep]
    file.remove(file.path(root.path,"Theta_opt",theta_to_delete))
    
    # Delete unwanted files in Evaluation
    eval_files <- list.files(file.path(root.path, "Evaluation"))
    eval_to_delete <- eval_files[basename(eval_files) != theta_keep]
    file.remove(file.path(root.path,"Evaluation", eval_to_delete))
    
    theta_final <- readRDS(file = file.path(root.path, "Theta_opt", paste0(beta, "_", lambda, ".rds"))) 
  }else{
    theta_final <- -sign(theta)*theta * 1e100
    res <- process_results(theta = theta_final, X =  X_test, A = A_test, 
                           Y = Y_test, Xi = Xi_test, mu0 = mu0_test, 
                           nu0 = nu0_test, prop_score = prop_score_test, 
                           lambda = lambda, alpha = alpha, beta = 0, 
                           centered = centered, conf_level = conf_level,
                           card_lambda = cardLambda, card_B = cardB)
    saveRDS(res, file=file.path(root.path,"Evaluation", 
                                paste0(0, "_", lambda,".rds")))
    attr(theta_final, "lambda") <- lambda
    attr(theta_final, "beta") <- 0
    saveRDS(theta_final, file = file.path(root.path, "Theta_opt", 
                                          paste0(0, "_", lambda, ".rds")))
  }
  psi_values <- make_psi(theta_final)(X_test)
  optimal_treatment_rule <- sigma_beta(psi_values, beta = attr(theta_final, "beta"))
    
  # Calculate proportion over 0.5
  prop_over_0.50 <- mean(optimal_treatment_rule > 0.5)
  if(prop_over_0.50<0.1){
    warning(sprintf(
      paste(
          "Only %.1f%% of the test set has an optimal treatment probability above 0.5.",
          "This may indicate that your tolerance for adverse events (alpha) is too strict.",
          "Consider relaxing it if treatment is being under-assigned."), 
        100 * prop_over_0.50))
    }
  return(theta_final)
}