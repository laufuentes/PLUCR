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
#' @param Jfold Number of folds for the main algorithm, needs to be set to 3L.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param SL.library Vector of libraries for training Super Learner (c("SL.mean","SL.glm","SL.ranger","SL.grf") by default).
#' @param tol A numeric scalar used as an early stopping criterion based on the RMSE between consecutive solutions (0.025 by default).
#' @param max_iter A numeric scalar specifying the maximum number of iterations (5 by default).
#' @param root.path Path to the folder where all results are to be saved.
#' @export
main_algorithm <- function(X, A, Y, Xi, 
                           Lambdas=seq(1, 8, by=1), alpha=0.1, precision=0.05,
                           B=c(0.05, 0.1, 0.25, 0.5), centered=FALSE,
                           Jfold=3, V=2L, SL.library=c("SL.mean","SL.glm","SL.ranger","SL.grf"), 
                           tol=0.025, max_iter=5, root.path){
  # Check whether the root.path exists and contains proper folder to save data
  if (!dir.exists(root.path)) {
    warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
    dir.create(root.path, recursive = TRUE)
  }
  # Subdirectories to check
  subdirs <- c("Mu.hat", "Nu.hat", "PS.hat", "Folds", "Images", "Intermediate", "Evaluation", "Theta_opt")
  
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
                                 cvControl = SuperLearner::SuperLearner.CV.control(V = JFold, shuffle = TRUE))
  saveRDS(folds, file=file.path(root.path,"Folds","folds.rds")) #Save primary outcome model
  # Check data
  checks<- PLUCR::check_data(Y, Xi, A, X, folds) 
  ########################################
  # Step 1: Compute nuisance parameters #
  ########################################
  # Primary outcome model 
  mu.hat.nj <- PLUCR::estimate_mu(Y=Y, A=A, X=X, folds=folds, SL.library = SL.library, V=2L)
  saveRDS(mu.hat.nj, file=file.path(root.path,"Mu.hat","mu.hat.nj.rds")) #Save primary outcome model
  # Adverse event model
  nu.hat.nj <- PLUCR::estimate_nu(Xi=Xi, A=A, X=X, folds=folds, SL.library = SL.library,V=2L)
  saveRDS(nu.hat.nj, file=file.path(root.path,"Nu.hat","nu.hat.nj.rds")) #Save adverse events outcome model
  # Propensity score
  ps.hat.nj <- PLUCR::estimate_ps(A=A, X=X, folds=folds, SL.library = SL.library,V=2L)
  saveRDS(ps.hat.nj, file=file.path(root.path,"PS.hat","ps.hat.nj.rds")) #Save propensity score outcome model
  
  ### Training data & functions
  # data
  X_train <- X[folds[[2]],]
  A_train <- A[folds[[2]]]
  Y_train <- Y[folds[[2]]]
  Xi_train <- Xi[folds[[2]]]
  # functions
  mu0_train <- function(a,x){mu.hat.nj(a=a,x=x,ff=1)}
  nu0_train <- function(a,x){nu.hat.nj(a=a,x=x,ff=1)}
  prop_score_train <- function(a,x){ps.hat.nj(a,x,1)}
  
  ###  Testing data & functions
  # data
  X_test <- X[folds[[3]],]
  A_test <- A[folds[[3]]]
  Y_test <- Y[folds[[3]]]
  Xi_test <- Xi[folds[[3]]]
  # functions
  mu0_test <- function(a,x){mu.hat.nj(a,x,3)}
  nu0_test <- function(a,x){nu.hat.nj(a,x,3)}
  prop_score_test <- function(a,x){ps.hat.nj(a,x,3)}
  
  ##########################################
  # Step 2: Training process: grid search  #
  ##########################################
  ##### 1.1- Start by testing \code{lambda}=0 
  saved <- FALSE
  combinations <- NULL
  lambda <- 0
  min_constraint_lambda0 <- Inf
  max_policy_value <- -Inf
  beta_0 <- NULL
  out <- PLUCR::Optimization_Estimation(mu0=mu0_train, nu0=nu0_train, prop_score=prop_score_train, 
                                        X=X_train, A=A_train, Y=Y_train, Xi=Xi_train, 
                                        lambda=0, alpha=alpha, precision=precision, beta=0.05, centered=centered, 
                                        tol=tol, max_iter=max_iter) #root.path=file.path(root.path,"Intermediate",paste0(0.05,"_",0))
  
  ##### 1.2- Check whether not considering your constraint satisfies already your constraint  
  theta_0 <- out$theta_collection[[length(out$theta_collection)]]
  attr(theta_0, "lambda") <- 0 
  for (beta in B){
    res_0 <- process_results(theta_0, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda=0, alpha,  beta, centered)
    saveRDS(out, file=file.path(root.path,"Intermediate",paste0(beta,"_",0,".rds")))
    saveRDS(res_0[[1]], file=file.path(root.path,"Evaluation",paste0(beta,"_",0,".rds")))
    
    # Loop to check constraint satisfaction
    if (res_0[[1]]$constraint < 0) {
      min_constraint_lambda0 <- res_0[[1]]$constraint
      
      if (res_0[[1]]$lwr_bound_policy_value > max_policy_value) {
        beta_0 <- beta
        max_policy_value <- res_0[[1]]$lwr_bound_policy_value
        saveRDS(theta_0, file = file.path(root.path, "Theta_opt", paste0(beta, "_", 0, ".rds")))
        saveRDS(res_0[[1]], file=file.path(root.path,"Evaluation",paste0(iteration,"_",beta,"_",0,".rds")))
        combinations <- rbind(combinations, c(beta_0, 0))
        saved <- TRUE}
    }else{
      if(res_0[[1]]$constraint<min_constraint_lambda0){
        min_constraint_lambda0 <- res_0[[1]]$constraint
        beta_0 <- beta
      }
    }
    ##### If your constraint was already satified with lambda=0 return
    if(res_0[[2]]<0){
      #/!\warning()
      return(theta_0)
    }
  }
  attr(theta_0, "beta") <- beta_0
  ##### If your constraint was not satified with lambda=0 goes onto step 2
  ##### 2.1- Test different lambda and beta combinations and save the optimal solutions satisfying the constraint 
    ### Training 
    for (beta in B){
      saved <- FALSE
      for (lambda in Lambdas){
        # Policy optimization
        out <- PLUCR::Optimization_Estimation(mu0=mu0_train, nu0=nu0_train, prop_score=prop_score_train, 
                                              X=X_train, A=A_train, Y=Y_train, Xi=Xi_train, 
                                              lambda=lambda, alpha=alpha, precision=precision, beta=beta, centered=centered, 
                                              tol=tol, max_iter=max_iter) # file.path(root.path,"Intermediate",paste0(beta,"_",lambda))
        theta_opt <- out$theta_collection[[length(out$theta_collection)]]
        ### Evaluating 
        res <- process_results(theta_opt, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda, alpha,  beta, centered)
        if (!saved && res[[1]]$constraint < 0) {
          saveRDS(res[[1]], file=file.path(root.path,"Evaluation", paste0(beta, "_", lambda,".rds")))
          saveRDS(out, file=file.path(root.path,"Intermediate",paste0(beta,"_",lambda,".rds")))
          saveRDS(theta_opt, file = file.path(root.path, "Theta_opt", paste0(beta, "_", lambda, ".rds")))
          saved <- TRUE
          combinations <- rbind(combinations, c(beta, lambda))
        }
        if(res[[2]]<0){
          break
        }
      }
    }
    optimal_combination <- get_opt_beta_lambda(combinations,root.path)
    theta_final <- readRDS(file = file.path(root.path, "Theta_opt", paste0(optimal_combination$beta, "_", optimal_combination$lambda, ".rds"))) 
  
    attr(theta_final, "lambda") <- optimal_combination$lambda 
    attr(theta_final, "beta") <- optimal_combination$beta
  return(list(theta_0,theta_final))
}