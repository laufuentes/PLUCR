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
#' @param Betas A vector of non-negative scalars controlling the sharpness of the treatment probability function (c(0.05, 0.5, 1, 2) by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param Jfold Number of folds for the main algorithm, needs to be set to 3L.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param SL.library Vector of libraries for training Super Learner (c("SL.mean","SL.glm","SL.ranger","SL.grf") by default).
#' @param root.path Path to the folder where all results are to be saved.
#' @export
main_algorithm <- function(X, A, Y, Xi, 
                           Lambdas=seq(1,8,by=1), alpha=0.1, precision=0.05,
                           Betas=c(0.05, 0.5, 1, 2), centered=FALSE,
                           Jfold=3, V=2L, 
                           SL.library=c("SL.mean","SL.glm","SL.ranger","SL.grf"), 
                           root.path){
  # Check whether the root.path exists and contains proper folder to save data
  if (!dir.exists(root.path)) {
    warning(sprintf("The directory '%s' does not exist. Creatung it...", root.path))
    dir.create(root.path, recursive = TRUE)
  }
  # Subdirectories to check
  subdirs <- c("Mu.hat", "Nu.hat", "PS.hat","Intermediate", "Evaluation", "Theta_opt")
  
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
  # Check data
  checks<- PLUCR::check_data(Y,Xi,A,X,folds) # /!\ changer ici s par fold 
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
  beta <-0.05 
  lambda <- 0
  out <- PLUCR::Optimization_Estimation(mu0_train, nu0_train, prop_score_train, 
                                        X_train, A_train, Y_train, Xi_train, 
                                        lambda=0, alpha, precision, beta, centered, 
                                        file.path(root.path,"Intermediate",paste0(0,"_",beta)))
  
  ##### 1.2- Check whether not considering your constraint satisfies already your constraint  
  theta_0 <- out$theta_collection[[length(out$theta_collection)]]
  res_0 <- process_results(theta_0, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda=0, alpha,  beta, centered)
  saveRDS(res_0[[1]], file=file.path(root.path,"Evaluation",paste0(0,"_",beta,".rds")))
  
  if (!saved && res_0[[1]]$constraint < 0) {
    saveRDS(theta_0, file = file.path(root.path, "Theta_opt", paste0(0, "_", beta, ".rds")))
    psi_0 <- make_psi(theta_opt)
    sigma_beta_0 <- sigma_beta(psi_0(X),beta=0.05,centered = FALSE)
    combinations <- rbind(combinations, c(beta, lambda))
    saved <- TRUE
  }
  ##### If your constraint was already satified with lambda=0 return
  if(res_0[[2]]<0){
    return(theta_0)
    }
  ##### If your constraint was not satified with lambda=0 goes onto step 2
  ##### 2.1- Test different lambda and beta combinations and save the optimal solutions satisfying the constraint 
  else{
    ### Training 
    for (beta in Betas){
      saved <- FALSE
      for (lambda in Lambdas){
        # Policy optimization
        out <- PLUCR::Optimization_Estimation(mu0_train, nu0_train, prop_score_train, 
                                              X_train, A_train, Y_train, Xi_train, 
                                              lambda, alpha, precision, beta, centered, file.path(root.path,"Intermediate",paste0(beta,"_",lambda)))
        theta_opt <- out$theta_collection[[length(out$theta_collection)]]
        ### Evaluating 
        res <- process_results(theta_opt, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda, alpha,  beta, centered)
        saveRDS(res[[1]], file=file.path(root.path,"Evaluation", paste0(beta, "_", lambda,".rds")))
        if (!saved && res[[1]]$constraint < 0) {
          saveRDS(theta_opt, file = file.path(root.path, "Theta_opt", paste0(beta, "_", lambda, ".rds")))
          saved <- TRUE
          psi_opt <- make_psi(theta_opt)
          sigma_beta_opt <- sigma_beta(psi_opt(X),beta=0.05,centered = FALSE)
          combinations <- rbind(combinations, c(beta, lambda))
        }
        if(res[[2]]<0){
          break
        }
      }
    }
    optimal_combination <- get_opt_beta_lambda(combinations,root.path)
    theta_final <- readRDS(file = file.path(root.path, "Theta_opt", paste0(optimal_combination$beta, "_", optimal_combination$lambda, ".rds"))) 
  }
  return(list(theta_0,theta_final))
}