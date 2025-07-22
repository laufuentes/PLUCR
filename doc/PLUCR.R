## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,        # show code
  message = TRUE,     # show messages from code
  warning = TRUE      # show warnings from code
)

## ----setup--------------------------------------------------------------------
library(PLUCR)

## ----eval=FALSE---------------------------------------------------------------
# devtools::install("path/to/PLUCR_package")

## ----eval=FALSE---------------------------------------------------------------
# library(PLUCR)
# library(SuperLearner)

## ----eval=FALSE---------------------------------------------------------------
# 
# n <- 2000  # Sample size
# exp <- PLUCR::data_gen(n, seed=2025) # changer scenario score de propension expit(-0.5*X3 + 0.5*1(X1>0.3))
# 
# # # Complete data (not available in real settings)
# df_complete <- exp[[1]]
# 
# # # Observational data (to be used for training)
# df_obs <- exp[[2]]

## ----eval=FALSE---------------------------------------------------------------
# # Covariates matrix
# X <- df_obs %>%
#   dplyr::select(starts_with("X."))%>%
#   as.matrix()
# 
# # Binary treatment assignment (0 = control, 1 = treated)
# A <- df_obs$Treatment
# 
# # Primary outcome, assumed to be scaled to [0, 1]
# Y <- df_obs$Y
# 
# # Binary indicator of adverse events (0=No adverse event, 1=adverse event)
# Xi <- df_obs$Xi

## ----eval=FALSE---------------------------------------------------------------
# # SuperLearner library used for estimating nuisance functions
# # You can customize this list. The default includes:
# # - SL.mean: simple mean predictor
# # - SL.glm: generalized linear model
# # - SL.ranger: random forest
# # - SL.grf: generalized random forest
# # - SL.xgboost: gradient boosting
# 
# SL.library <- c("SL.mean", "SL.glm", "SL.ranger", "SL.grf", "SL.xgboost")
# 
# # Root directory where output files and intermediate results will be saved
# root.path <- "results"
# 
# # Constraint tolerance (alpha): maximum allowed increase in adverse event risk
# # Example: 0.15 means the policy may increase adverse event probability by up to 15%
# alpha.tolerance <- 0.15
# 
# # Optimization precision: controls granularity of the policy search
# # Smaller values lead to finer approximations but require more computation
# optimization.precision <- 0.05

## ----eval=FALSE---------------------------------------------------------------
# output <- PLUCR::main_algorithm(X=X, A=A, Y=Y, Xi=Xi,
#                                 alpha=alpha.tolerance, precision=optimization.precision,
#                                 SL.library= SL.library, root.path=root.path)

## ----eval=FALSE---------------------------------------------------------------
# n <- nrow(X)
# 
# # Folds for cross-validation
# folds <- SuperLearner::CVFolds(n,
#                                id = NULL,
#                                Y = Y,
#                                cvControl = SuperLearner::SuperLearner.CV.control(V = JFold,
#                                                                                  shuffle =TRUE))
# # Check data
# checks<- PLUCR::check_data(Y, Xi, A, X, folds)

## ----eval=FALSE---------------------------------------------------------------
# # Estimate primary outcome model E[Y | A, X]
# mu.hat.nj <- PLUCR::estimate_mu(
#   Y = Y, A = A, X = X, folds = folds,
#   SL.library = SL.library, V = 2L
# )
# saveRDS(mu.hat.nj, file = file.path(root.path, "Mu.hat", "mu.hat.nj.rds"))
# 
# # Estimate adverse event model E[Xi | A, X]
# nu.hat.nj <- PLUCR::estimate_nu(
#   Xi = Xi, A = A, X = X, folds = folds,
#   SL.library = SL.library, V = 2L
# )
# saveRDS(nu.hat.nj, file = file.path(root.path, "Nu.hat", "nu.hat.nj.rds"))
# 
# # Estimate propensity score model P[A | X]
# ps.hat.nj <- PLUCR::estimate_ps(
#   A = A, X = X, folds = folds,
#   SL.library = SL.library, V = 2L
# )
# saveRDS(ps.hat.nj, file = file.path(root.path, "PS.hat", "ps.hat.nj.rds"))

## ----eval=FALSE---------------------------------------------------------------
# # Training functions & data (first and second folds)
# # functions (first fold)
# mu0_train <- function(a,x){mu.hat.nj(a=a,x=x,ff=1)}
# nu0_train <- function(a,x){nu.hat.nj(a=a,x=x,ff=1)}
# prop_score_train <- function(a,x){ps.hat.nj(a,x,1)}
# 
# # data (second fold)
# X_train <- X[folds[[2]],]
# A_train <- A[folds[[2]]]
# Y_train <- Y[folds[[2]]]
# Xi_train <- Xi[folds[[2]]]
# 
# 
# # Testing functions & data (third fold)
# # functions
# mu0_test <- function(a,x){mu.hat.nj(a,x,3)}
# nu0_test <- function(a,x){nu.hat.nj(a,x,3)}
# prop_score_test <- function(a,x){ps.hat.nj(a,x,3)}
# 
# # data (third fold)
# X_test <- X[folds[[3]],]
# A_test <- A[folds[[3]]]
# Y_test <- Y[folds[[3]]]
# Xi_test <- Xi[folds[[3]]]

## ----eval=FALSE---------------------------------------------------------------
# # Start by testing the unconstrained case with lambda = 0
# saved <- FALSE  # Flag to track whether a valid result has already been saved
# combinations <- NULL  # Matrix to store (beta, lambda) pairs that satisfy the constraint
# beta <- 0.05  # Arbitrary beta value since constraint is not evaluated yet
# lambda <- 0
# 
# # Run optimization with lambda = 0 (no constraint)
# out <- PLUCR::Optimization_Estimation(mu0_train, nu0_train, prop_score_train,
#                                         X_train, A_train, Y_train, Xi_train,
#                                         lambda=0, alpha, precision, beta, centered,
#                                         file.path(root.path,"Intermediate",paste0(0,"_",beta)))
# 
# # Extract the final policy from the iteration
# theta_0 <- out$theta_collection[[length(out$theta_collection)]]
# # Evaluate the policy on the test data
# res_0 <- process_results(theta_0,
#                          X_test, A_test, Y_test, Xi_test,
#                          mu0_test, nu0_test, prop_score_test,
#                          lambda=0, alpha,  beta, centered)
# # Save the evaluation results
# saveRDS(res_0[[1]], file=file.path(root.path,"Evaluation",paste0(0,"_",beta,".rds")))
# 
# # If the constraint is satisfied, save the corresponding policy and mark as saved
# if (!saved && res_0[[1]]$constraint < 0) {
#     saveRDS(theta_0, file = file.path(root.path, "Theta_opt", paste0(0, "_", beta, ".rds")))
#     psi_0 <- make_psi(theta_opt)
#     sigma_beta_0 <- sigma_beta(psi_0(X), beta=0.05, centered = FALSE)
#     combinations <- rbind(combinations, c(beta, lambda))
#     saved <- TRUE
# }
# # If the upper bound of the constraint's CI is below zero, stop here
# if(res_0[[2]]<0){
#     print(theta_0)
# }else{
#     # Begin training over the grid
#     ##for (beta in Betas){
#       for (lambda in Lambdas){
#        # Run alternated procedure for each beta and lambda combination
#         out <- PLUCR::Optimization_Estimation(mu0_train, nu0_train, prop_score_train,
#                                               X_train, A_train, Y_train, Xi_train,
#                                               lambda, alpha, precision, beta, centered,
#                                             file.path(root.path,"Intermediate",paste0(lambda,"_",beta)))
#         # Extract final policy
#         theta_opt <- out$theta_collection[[length(out$theta_collection)]]
#          # Evaluate the policy on test set
#         res <- process_results(theta_opt,
#                                X_test, A_test, Y_test, Xi_test,
#                                mu0_test, nu0_test, prop_score_test,
#                                lambda, alpha, beta, centered)
#         # Save evaluation results
#         saveRDS(res[[1]], file=file.path(root.path,"Evaluation", paste0(lambda,"_",beta,".rds")))
# 
#         # If constraint is satisfied for the first time, save the optimal policy
#         if (!saved && res[[1]]$constraint < 0) {
#           saveRDS(theta_opt,
#                   file = file.path(root.path, "Theta_opt", paste0(lambda, "_", beta, ".rds")))
#           saved <- TRUE
#           psi_opt <- make_psi(theta_opt)
#           sigma_beta_opt <- sigma_beta(psi_opt(X), beta=0.05, centered = FALSE)
#           combinations <- rbind(combinations, c(beta, lambda))
#         }
#         # Stop search if constraint's upper CI bound is negative
#         if(res[[2]]<0){
#           break
#         }
# 
#       #}
#       }
#     theta_final <- readRDS(file = file.path(root.path, "Theta_opt",
#                                             paste0(combinations[,2], "_", combinations[,1], ".rds")))
#   }
# print(list(theta_0, theta_final))

