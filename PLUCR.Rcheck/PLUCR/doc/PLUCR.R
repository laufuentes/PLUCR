## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,        # show code
  message = TRUE,     # show messages from code
  warning = TRUE      # show warnings from code
)

## ----eval=FALSE---------------------------------------------------------------
# devtools::install_local(
#     "path/to/downloaded_file.tar.gz",
#     dependencies = TRUE # Install dependencies from CRAN
# )

## ----eval=FALSE---------------------------------------------------------------
# library(PLUCR)
# library(SuperLearner)
# library(grf)
# library(dplyr)
# library(tidyr)
# library(ggplot2)

## ----eval=FALSE---------------------------------------------------------------
# n <- 3000  # Sample size
# ncov <- 10L # Number of covariates
# scenario_mu <- "Linear"
# scenario_nu <- "Linear"
# 
# # Generate synthetic scenario
# exp <- PLUCR::generate_data(n,
#                             ncov,
#                             scenario_mu = scenario_mu,
#                             scenario_nu = scenario_nu,
#                             seed=2025)
# 
# # # Observational data (to be used for training)
# df_obs <- exp[[2]]

## ----eval=FALSE---------------------------------------------------------------
# # Extract normalized numerical covariates as a matrix
# # ---------------------------------------------------------------
# # Reminder:
# # - Ensure categorical variables are one-hot encoded
# #   (dummy variables) prior to inclusion.
# # ---------------------------------------------------------------
# 
# X <- df_obs %>%
#   dplyr::select(starts_with("X."))%>%
#   as.matrix()
# 
# # ---------------------------------------------------------------
# # Optional: Normalize covariates
# # Applying phi(X) rescales each column of the covariate matrix X
# # (e.g., centering and scaling) to improve numerical stability and
# # ensure comparability across features.
# # ---------------------------------------------------------------
# # X <- phi(X)
# 
# # Binary treatment assignment (0 = control, 1 = treated)
# A <- df_obs$Treatment
# 
# # Primary outcome, assumed to be scaled to [0, 1]
# Y <- df_obs$Y
# 
# # Note: For a normalization of the primary outcome, you can use:
# # Y <- (df_obs$Y - min(df_obs$Y))/(max(df_obs$Y)-min(df_obs$Y))
# 
# # Binary indicator of adverse events
# # (0 = No adverse event, 1 = adverse event)
# Xi <- df_obs$Xi

## ----eval=FALSE---------------------------------------------------------------
# # SuperLearner library used for estimating nuisance functions
# # You can customize this list. The default includes:
# # - SL.mean: simple mean predictor
# # - SL.glm: generalized linear model
# # - SL.ranger: random forest
# # - SL.grf: generalized random forest
# # - SL.xgboost: gradient boosting
# SL.library <- c("SL.mean", "SL.glm", "SL.ranger", "SL.grf", "SL.xgboost")
# 
# # Insert a root directory where output files, results and images will be saved
# root.path <- "results"
# 
# # An increasing sequence of non-negative numeric scalars controlling
# # the penalty for violating the constraint (seq(1,8, by=1) by default).
# Lambda_candidates <- seq(1, 10, by=1)
# 
# # A vector of non-negative scalars controlling the sharpness of the
# # treatment probability function (c(0, 0.05, 0.1, 0.25, 0.5) by default)
# Beta_candidates <- c(0, 0.05, 0.1, 0.25, 0.5)
# 
# # Constraint tolerance (alpha): maximum allowed increase in adverse event risk
# # Example: 0.2 means the policy may increase adverse event probability by up to 20%
# alpha.tolerance <- 0.2
# 
# # Optimization precision: Controls the granularity of the policy search.
# # Smaller values lead to finer approximations but require more computation (default is 0.025).
# # For a quick toy example, you may prefer using 0.05 instead.
# optimization.precision <- 0.025
# 
# # Parameter indicating whether to center the treatment policy to ensure that
# # probability of treatment assignment for null treatment effect is 0.5.
# centered <- FALSE

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Run the Main Algorithm to Learn the Optimal Treatment Rule
# # ===============================================================
# 
# output <- PLUCR::main_algorithm(X=X, A=A, Y=Y, Xi=Xi,
#                                 Lambdas = Lambda_candidates,
#                                 alpha=alpha.tolerance,
#                                 B= Beta_candidates,
#                                 precision=optimization.precision,
#                                 centered=centered,
#                                 SL.library = SL.library,
#                                 root.path=root.path)
# 
# # ---------------------------------------------------------------
# # Extract learned parameters from output
# # Depending on configuration, the algorithm may return:
# #   - A list with (1) initial estimate theta_0 and (2) final estimate theta_final
# #   - A single theta object if constraint satisfied at lambda=0
# # ---------------------------------------------------------------
# 
# if(is.list(output)){
#   theta_0 <- output[[1]]
#   theta_final <- output[[2]]
# }else{
#   theta_final <- output
# }

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Derive the Optimal Treatment Rule from Final Parameters
# # ===============================================================
# # ---------------------------------------------------------------
# # Step 1. Construct the decision function `psi` using the `make_psi`
# # It maps covariates X to [-1,1] scores for treatment assignment.
# # ---------------------------------------------------------------
# psi_final <- PLUCR::make_psi(theta_final)
# 
# # Obtain the optimal treatment rule by applying `sigma_beta` to the optimal `beta` associated to `theta_final`
# # ---------------------------------------------------------------
# # Step 2. Translate `psi(X)` into treatment probabilities
# # Apply the treatment probability function `sigma_beta`,
# # controlled by beta.
# # ---------------------------------------------------------------
# optimal_treatment_probability <- PLUCR::sigma_beta(psi_final(X), beta=attr(theta_final, "beta"))

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Evaluate the Optimal Treatment Rule using normalized data
# # ===============================================================
# # ---------------------------------------------------------------
# # Step 1. Load pre-required objects (algorithm outputs)
# # ---------------------------------------------------------------
# folds <- readRDS(file.path(root.path,"Folds","folds.rds")) # Folds for splitting data
# mu.hat.nj <- readRDS(file.path(root.path,"Mu.hat","mu.hat.nj.rds")) # Primary outcome model
# nu.hat.nj <- readRDS(file.path(root.path,"Nu.hat","nu.hat.nj.rds")) # Adverse event model
# ps.hat.nj <- readRDS(file.path(root.path,"PS.hat","ps.hat.nj.rds")) # Propensity score model
# 
# # ---------------------------------------------------------------
# # Step 2. Extract the test set for policy evaluation
# # (thereby ensuring unbiased estimates of treatment performance)
# # ---------------------------------------------------------------
# X_test <- X[folds[[3]],]
# A_test <- A[folds[[3]]]
# Y_test <- Y[folds[[3]]]
# Xi_test <- Xi[folds[[3]]]
# 
# # ---------------------------------------------------------------
# # Step 3. Extract test nuisance functions
# # ---------------------------------------------------------------
# mu0_test <- function(a,x){mu.hat.nj(a,x,3)}
# nu0_test <- function(a,x){nu.hat.nj(a,x,3)}
# prop_score_test <- function(a,x){ps.hat.nj(a,x,3)}
# 
# # ---------------------------------------------------------------
# # Step 4. Evaluate the learned optimal treatment rule
# # ---------------------------------------------------------------
# eval_plucr <- process_results(theta=theta_final, X=X_test, A=A_test, Y=Y_test,
#                               Xi=Xi_test, mu0=mu0_test, nu0=nu0_test,
#                               prop_score=prop_score_test,
#                               lambda=attr(theta_final, "lambda"),
#                               alpha=alpha.tolerance,
#                               beta=attr(theta_final, "beta"),
#                               centered=centered)
# 
# # ---------------------------------------------------------------
# # Optional: Oracular evaluation (synthetic scenarios only)
# # Uses the known data-generating process.
# # ---------------------------------------------------------------
# # Extract oracular features (only available in synthetic settings)
# df_complete <- exp[[1]] # Complete data including counterfactuals c(Y(0),Y(1), xi(0), xi(1))
# delta_mu <- exp[[3]] # Oracular delta_mu function (treatment effect over Y conditioned on X)
# delta_nu <- exp[[4]] # Oracular delta_nu function (treatment effect over xi conditioned on X)
# 
# # Oracular evaluation
# eval_plucr_oracular <- PLUCR::oracular_process_results(theta=theta_final, ncov=ncov,
#                                               scenario_mu = scenario_mu,
#                                               scenario_nu = scenario_nu,
#                                               lambda = attr(theta_final, "lambda"),
#                                               alpha = alpha.tolerance,
#                                               beta= attr(theta_final, "beta"))

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Learn Decision Threshold and Derive Treatment Assignments
# # ===============================================================
# # ---------------------------------------------------------------
# # Step 1. Learn the decision threshold t
# # ---------------------------------------------------------------
# t <- learn_threshold(theta_final, X, A, Y, Xi, folds, alpha.tolerance)
# 
# # ---------------------------------------------------------------
# # Step 2. Convert treatment probabilities into binary decisions
# # Assign treatment if the individualized probability exceeds t.
# # ---------------------------------------------------------------
# optimal_treatment_decision <- ifelse(optimal_treatment_probability>t, 1, 0)

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Evaluate the Optimal Treatment Rule using non-normalized data
# # ===============================================================
# # ---------------------------------------------------------------
# # Step 1. Train pre-required objects (algorithm outputs)
# # ---------------------------------------------------------------
# Y_non_normalized <- # Insert non scaled primary outcome
# covariates <- # Insert X prior to `phi(X)`
# mu.hat.nj <- estimate_real_valued_mu(Y= Y_non_normalized, A= A, X=covariates, folds, SL.library = SL.library, V=2L)
# nu.hat.nj <- estimate_nu(Xi= Xi, A= A, X=covariates, folds, SL.library = SL.library, V=2L, threshold = 0.01)
# 
# # Nuisances trained on testing subset for unbiased evaluation
# mu0_test <- function(a,x){mu.hat.nj(a,x,ff=3)}
# nu0_test <- function(a,x){nu.hat.nj(a,x,ff=3)}
# 
# # Data testing subset
# covariates_test <- covariates[folds[[3]],]
# Y_non_test <- Y_non_normalized[folds[[3]]]
# proba_pi <- optimal_treatment_probability[folds[[3]]]
# pi_decision <- optimal_treatment_decision[folds[[3]]]
# 
# # Step 2. Evaluate probabilistic policy and derived decision
# # ---------------------------------------------------------------
# # Evaluation of optimal treatment probabilities
# policy_value_proba <- V_Pn(proba_pi,
#                            y1 = mu0_test(a = rep(1, nrow(covariates_test)), covariates_test),
#                            y0 = mu0_test(a = rep(0, nrow(covariates_test)), covariates_test))
# 
# constraint_value_proba <- mean(proba_pi*(nu0_test(a = rep(1, nrow(covariates_test)), covariates_test) - nu0_test(a = rep(0, nrow(covariates_test)), covariates_test))) - alpha.tolerance
# 
# print(c(policy_value_proba, constraint_value_proba))
# 
# # Evaluation of optimal treatment decisions derived from probabilities
# policy_value_decision <- V_Pn(pi_decision,
#                               y1 = mu0_test(a = rep(1, nrow(covariates_test)), covariates_test),
#                               y0 = mu0_test(a = rep(0, nrow(covariates_test)), covariates_test))
# 
# constraint_value_decision <- mean(pi_decision*(nu0_test(a = rep(1, nrow(covariates_test)), covariates_test) - nu0_test(a = rep(0, nrow(covariates_test)), covariates_test))) - alpha.tolerance
# 
# print(c(policy_value_decision, constraint_value_decision))

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Run the Naive Algorithm to Learn the Optimal Treatment Rule
# # ===============================================================
# naive_output <- naive_approach_algorithm(X=X, A=A, Y=Y, Xi=Xi,
#                                          folds= folds,
#                                          mu0_train = mu0_train,
#                                          mu0_test = mu0_test,
#                                          nu0_train = nu0_train,
#                                          nu0_test = nu0_test,
#                                          prop_score_train = prop_score_train,
#                                          prop_score_test = prop_score_test,
#                                          alpha = alpha.tolerance,
#                                          centered=centered,
#                                          precision=optimization.precision,
#                                          root.path = root.path)
# 
# # ---------------------------------------------------------------
# # Extract learned parameters from output
# # Depending on configuration, the algorithm may return:
# #   - A list with (1) initial estimate theta_0 and (2) final estimate theta_final
# #   - A single theta object if constraint satisfied at lambda=0
# # ---------------------------------------------------------------
# if(is.list(naive_output)){
#  theta_naive <- naive_output[[2]]
# }else{
#   theta_naive <- naive_output[[1]]
# }
# 
# # =================================================================
# # Evaluate the Optimal Treatment Rule according to Naive Approach
# # =================================================================
# eval_naive <- process_results(theta=theta_naive, X=X_test, A=A_test, Y=Y_test,
#                               Xi=Xi_test, mu0=mu0_test, nu0=nu0_test,
#                               prop_score=prop_score_test,
#                               lambda=attr(theta_naive, "lambda"),
#                               alpha=alpha.tolerance,
#                               beta=attr(theta_naive, "beta"),
#                               centered=centered)
# 
# # ---------------------------------------------------------------
# # Optional: Oracular evaluation (synthetic scenarios only)
# # Uses the known data-generating process.
# # ---------------------------------------------------------------
# eval_naive_oracular <- PLUCR::oracular_process_results(theta=theta_naive, ncov=ncov,
#                                               scenario_mu = scenario_mu,
#                                               scenario_nu = scenario_nu,
#                                               lambda = attr(theta_naive, "lambda"),
#                                               alpha = alpha.tolerance,
#                                               beta= attr(theta_naive, "beta"))

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Run the Orcular Algorithm to Learn the Optimal Treatment Rule
# # ===============================================================
# oracular_output <- oracular_approach_algorithm(X=X, A=A, Y=Y, Xi=Xi,
#                                                folds=folds,  ncov=ncov,
#                                                delta_Mu=delta_mu, delta_Nu= delta_nu,
#                                                alpha = alpha.tolerance,
#                                                centered=centered,
#                                                precision = optimization.precision,
#                                                scenario_mu=scenario_mu,
#                                                scenario_nu=scenario_nu,
#                                                root.path = root.path)
# 
# # ---------------------------------------------------------------
# # Extract learned parameters from output
# # Depending on configuration, the algorithm may return:
# #   - A list with (1) initial estimate theta_0 and (2) final estimate theta_final
# #   - A single theta object if constraint satisfied at lambda=0
# # ---------------------------------------------------------------
# if(is.list(oracular_output)){
#  theta_oracular <- oracular_output[[2]]
# }else{
#   theta_oracular <- oracular_output[[1]]
# }
# 
# # =================================================================
# # Evaluate oracularly the Optimal Treatment Rule
# # =================================================================
# eval_oracular <- oracular_process_results(theta=theta_oracular, ncov=ncov,
#                                           scenario_mu = scenario_mu,
#                                           scenario_nu = scenario_nu,
#                                           lambda = attr(theta_oracular, "lambda"),
#                                           alpha = alpha.tolerance,
#                                           beta= attr(theta_oracular, "beta"))

## ----eval=FALSE---------------------------------------------------------------
# # ---------------------------------------------------------------
# # Plot treatment effect functions:
# #   - delta_Mu: for primary outcome function
# #   - delta_Nu: for adverse event function
# # The plot helps to interpret how treatment impacts outcomes
# # under the chosen scenario.
# # ---------------------------------------------------------------
# fig <- PLUCR::synthetic_data_plot(delta_Mu=delta_mu,
#                                   delta_Nu=delta_nu,
#                                   root.path=root.path,
#                                   name=paste0(scenario_mu,"-", scenario_nu))
# 
# print(fig)

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Visualize Learned Treatment Policies:
# # Plot individualized treatment probabilities as a function of
# # two selected covariates (X.1 and X.2) to enable 2-D visualization
# # of the treatment probabilities.
# # ===============================================================
# 
# if(is.list(output)){
#  PLUCR::visual_treatment_plot(make_psi(theta_0)(X),
#                               lambda = attr(theta_0, "lambda"),
#                               beta = attr(theta_0, "beta"),
#                               centered = centered,
#                               Var_X_axis = df_obs$X.1, Var_Y_axis = df_obs$X.2,
#                               root.path = root.path, name = "Initial")
# }
# 
#  PLUCR::visual_treatment_plot(make_psi(theta_naive)(X),
#                               lambda = attr(theta_naive, "lambda"),
#                               beta = attr(theta_naive, "beta"),
#                               centered = centered,
#                               Var_X_axis = df_obs$X.1, Var_Y_axis = df_obs$X.2,
#                               root.path = root.path, name = "Naive")
# 
#  PLUCR::visual_treatment_plot(make_psi(theta_final)(X),
#                               lambda = attr(theta_final, "lambda"),
#                               beta = attr(theta_final, "beta"),
#                               centered = centered,
#                               Var_X_axis = df_obs$X.1, Var_Y_axis = df_obs$X.2,
#                               root.path = root.path, name = "Final")
# 
#  PLUCR::visual_treatment_plot(make_psi(theta_oracular)(X),
#                               lambda = attr(theta_oracular, "lambda"),
#                               beta = attr(theta_oracular, "beta"),
#                               centered = centered,
#                               Var_X_axis = df_obs$X.1, Var_Y_axis = df_obs$X.2,
#                               root.path = root.path, name = "Oracular")

## ----eval=FALSE---------------------------------------------------------------
# # ===============================================================
# # Visualize metrics for different outcomes
# # ===============================================================
# eval_0 <- process_results(theta_0, X_test, A_test, Y_test, Xi_test, mu0_test, nu0_test, prop_score_test, lambda=attr(theta_0, "lambda"), alpha.tolerance,  beta=attr(theta_0, "beta"), centered)[[1]]
# 
# if(is.list(output)){
#   data <- tibble(
#   method = c("Theta.0", "Theta.naive.PLUCR", "Theta.PLUCR", "Theta.Oracle.PLCUR"),
#   policy_value = c(eval_0$policy_value,
#                    eval_naive_oracular$policy_value,
#                    eval_plucr_oracular$policy_value,
#                    eval_oracular$policy_value),
#   constraint = c(eval_0$constraint,
#                  eval_naive_oracular$constraint,
#                  eval_plucr_oracular$constraint,
#                  eval_oracular$constraint))
# }else{
#   data <- tibble(
#   method = c("Theta_naive", "Theta_final", "Theta_oracular"),
#   policy_value = c(eval_naive_oracular$policy_value,
#                    eval_plucr_oracular$policy_value,
#                    eval_oracular$policy_value),
#   constraint = c(eval_naive_oracular$constraint,
#                  eval_plucr_oracular$constraint,
#                  eval_oracular$constraint))
# }
# 
# # Call the function
# plot_metric_comparison(data, metrics = select(data, -"method") %>% colnames(),
#                        techniques = data$method, root.path = root.path)

## ----eval=FALSE---------------------------------------------------------------
#  n <- nrow(X)
# 
# # Folds for cross-validation
# folds <- SuperLearner::CVFolds(n,
#                                id = NULL,
#                                Y = Y,
#                                cvControl = SuperLearner::SuperLearner.CV.control(V = Jfold,
#                                                                                  shuffle = TRUE))
# 
# saveRDS(folds, file=file.path(root.path,"Folds","folds.rds"))
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
# beta_0 <- NULL  # Arbitrary beta value since constraint is not evaluated yet
# lambda <- 0
# 
# # Run optimization with lambda = 0 (no constraint)
# out <- PLUCR::Optimization_Estimation(mu0=mu0_train, nu0=nu0_train,
#                                       prop_score=prop_score_train,
#                                       X=X_train, A=A_train, Y=Y_train, Xi=Xi_train,
#                                       lambda=0, alpha=alpha, precision=precision, beta=0.05,
#                                       centered=centered, tol=tol, max_iter=max_iter)
# 
# # Extract the final policy from the iteration
# theta_0 <- out$theta_collection[[length(out$theta_collection)]]
# attr(theta_0, "lambda") <- 0  # Save the lambda as an attribute
# 
# # Evaluate the policy on the test data for all possible betas
# for (beta in B){
#     res_0 <- process_results(theta_0, X_test, A_test, Y_test, Xi_test,
#                              mu0_test, nu0_test, prop_score_test, lambda=0,
#                              alpha,  beta, centered)
#     # Loop to check constraint satisfaction
#     if (res_0[[1]]$constraint < 0) {
#       min_constraint_lambda0 <- res_0[[1]]$constraint
# 
#       if (res_0[[1]]$lwr_bound_policy_value > max_policy_value) {
#         beta_0 <- beta
#         max_policy_value <- res_0[[1]]$lwr_bound_policy_value
#         saveRDS(theta_0, file = file.path(root.path, "Theta_opt", paste0(beta, "_", 0, ".rds")))
#          attr(theta_0, "lambda") <- lambda
#           attr(theta_0, "beta") <- beta
#         saveRDS(out, file=file.path(root.path,"Intermediate",paste0(beta,"_",0,".rds")))
#         saveRDS(res_0[[1]], file=file.path(root.path,"Evaluation",paste0(beta,"_",0,".rds")))
#         combinations <- rbind(combinations, c(beta_0, 0))
#         saved <- TRUE}
#     }else{
#       if(res_0[[1]]$constraint<min_constraint_lambda0){
#         min_constraint_lambda0 <- res_0[[1]]$constraint
#         beta_0 <- beta
#       }
#   }
#     if(res_0[[2]]<0){
#       warning(sprintf(paste("The constraint was already satisfied for lambda=0.")))
#       attr(theta_0, "beta") <- beta_0
#       return(theta_0)
#     }
#   }
#   attr(theta_0, "beta") <- beta_0 # Save the beta as an attribute
# 
#   # Begin training over the grid: find for each beta value the optimal lambda
#  for (beta in B){
#       saved <- FALSE
#       for (lambda in Lambdas){
#        # Run alternated procedure for each beta and lambda combination
#          out <- PLUCR::Optimization_Estimation(mu0=mu0_train, nu0=nu0_train,
#                                                prop_score=prop_score_train,
#                                                X=X_train, A=A_train, Y=Y_train, Xi=Xi_train,
#                                                lambda=lambda, alpha=alpha, precision=precision,
#                                                beta=beta, centered=centered,
#                                                tol=tol, max_iter=max_iter)
#         # Extract final policy
#         theta_opt <- out$theta_collection[[length(out$theta_collection)]]
#         ### Evaluating
#         res <- process_results(theta_opt, X_test, A_test, Y_test, Xi_test,
#                                mu0_test, nu0_test, prop_score_test, lambda,
#                                alpha,  beta, centered)
#         if (!saved && res[[1]]$constraint < 0) {
#           saveRDS(res[[1]], file=file.path(root.path,"Evaluation",
#                                            paste0(beta, "_", lambda,".rds")))
#           saveRDS(out, file=file.path(root.path,"Intermediate",
#                                       paste0(beta,"_",lambda,".rds")))
#           saveRDS(theta_opt, file = file.path(root.path, "Theta_opt",
#                                               paste0(beta, "_", lambda, ".rds")))
#           attr(theta_opt, "lambda") <- lambda
#           attr(theta_opt, "beta") <- beta
#           saved <- TRUE
#           combinations <- rbind(combinations, c(beta, lambda))
#         }
#         # Stop search if constraint's upper CI bound is negative
#         if(res[[2]]<0){
#           break
#         }
#       }
#  }
# # Select the optimal combination (beta, lambda)
# files_del <- file.path(root.path,"Intermediate")
# unlink(files_del, recursive = TRUE)
# 
# # Select the optimal combination (beta, lambda)
# optimal_combination <- get_opt_beta_lambda(combinations,root.path)
# beta <- optimal_combination$beta
# lambda <- optimal_combination$lambda
# theta_keep <- paste0(beta, "_", lambda, ".rds")
# 
# # Delete unwanted files in Theta_opt
# theta_files <- list.files(file.path(root.path, "Theta_opt"))
# theta_to_delete <- theta_files[basename(theta_files) != theta_keep]
# file.remove(file.path(root.path,"Theta_opt",theta_to_delete))
# 
# # Delete unwanted files in Evaluation
# eval_files <- list.files(file.path(root.path, "Evaluation"))
# eval_to_delete <- eval_files[basename(eval_files) != theta_keep]
# file.remove(file.path(root.path,"Evaluation", eval_to_delete))
# 
# # Load the corresponding theta for the optimal (beta, lambda) combination
# theta_final <- readRDS(file = file.path(root.path, "Theta_opt",
#                                             paste0(optimal_combination$beta, "_",
#                                                    optimal_combination$lambda, ".rds")))
# # Save the optimal combiation (lambda, beta) as an attribute
# attr(theta_final, "lambda") <- optimal_combination$lambda
# attr(theta_final, "beta") <- optimal_combination$beta
# 
# # Compute the optimal treatment probabilities with sigma_beta
# psi_values <- make_psi(theta_final)(X_test)
# optimal_treatment_rule <- sigma_beta(psi_values, beta = attr(theta_final, "beta"))
# 
# # Calculate proportion over 0.5 and warn the user if the treatment probabilities are too low.
# prop_over_0.50 <- mean(optimal_treatment_rule > 0.5)
# if(prop_over_0.50<0.1){
#   warning(sprintf(
#         paste("Only %.1f%% of the test set has an optimal treatment probability above 0.5.",
#           "This may indicate that your tolerance for adverse events (alpha) is too strict.",
#           "Consider relaxing it if treatment is being under-assigned."), 100 * prop_over_0.50))
# }
# 
# print(list(theta_0, theta_final))

