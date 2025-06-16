#' Estimate mu 
#'
#' This function trains conditional mean of primary outcome models for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in [0, 1]).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param folds A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}. 
#' @param SL.library Vector of libraries for training SuperLearner.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-3 by default). Constrains estimation to $[threshold, 1 - threshold]$.
#' @return A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X)
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' A <- stats::rbinom(100, 1, 0.5)
#' Y <- runif(100)
#' JFold <- 3
#' folds <- folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y, cvControl = SuperLearner::SuperLearner.CV.control(V = JFold, shuffle = TRUE))
#' SL.library <- c("SL.glm", "SL.mean", "SL.ranger")
#' mu_functions <- estimate_mu(Y, A, X, folds)
#' # Apply a function from the list to new data:
#' mu_functions(1,X[1:5, ],1)  # Predict treated outcomes
#' }
#' @export
estimate_mu <- function(Y, A, X, folds, SL.library, V = 2L, threshold = 1e-3) {
  if (!(threshold<5*1e-2 & threshold>0)) {
    msg_threshold <- paste("Threshold:", threshold, "is not small enough or not positive")
    warning(msg_threshold)
  }
  if (!(V<5 & V>1)) {
    msg_folds <- paste("Number of folds:",V, "is either smaller than 1 or greater than 10")
    warning(msg_folds)
  }
  if (missing(SL.library)) {
    SL.library <- c("SL.mean", "SL.glm")
  }
  if (!is.data.frame(X)){
    X <- as.data.frame(X)
  }
  cvControl <- SuperLearner::SuperLearner.CV.control(V = V)
  method <- "method.NNLS"
  n <- nrow(X)
  objects <- vector("list", length(folds))
  for (ff in 1:length(folds)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(Y[train_indices],
                                                data.frame(X, A = A)[train_indices, , drop = FALSE],
                                                SL.library = SL.library,
                                                family = gaussian(),
                                                cvControl = cvControl, method = method)
  }
  mu <- function(a, x, ff) {
    if (length(a) != nrow(x)) stop("Length of 'a' must match number of rows in 'x'")
    out <- rep(NA, length(a))
    out <- stats::predict(objects[[ff]],newdata = data.frame(x,A=a))$pred
    out <- pmax(threshold, pmin(1 - threshold, out))
    return(out)
  }
  
  return(mu)
}

#' Estimate nu 
#'
#' This function trains conditional mean of adverse event outcome models for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param Xi A numeric vector or matrix of adverse events outcomes.
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param folds A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}. 
#' @param SL.library Vector of libraries for training SuperLearner.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-3 by default). Constrains estimation to $[threshold, 1 - threshold]$.
#' @return A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X)
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' A <- stats::rbinom(100, 1, 0.5)
#' Xi <- rbinom(100, n=1, p=0.5)
#' JFold <- 3
#' folds <- folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y, cvControl = SuperLearner::SuperLearner.CV.control(V = JFold, shuffle = TRUE))
#' SL.library <- c("SL.glm", "SL.mean", "SL.ranger")
#' nu_functions <- estimate_nu(Xi, A, X, folds)
#' # Apply a function from the list to new data:
#' nu_functions(1,X[1:5, ],1)  # Predict treated outcomes
#' }
#' @export
estimate_nu <- function(Xi, A, X, folds, SL.library, V = 2L, threshold = 1e-3) {
  if (!(threshold<5*1e-2 & threshold>0)) {
    msg_threshold <- paste("Threshold:", threshold, "is not small enough or not positive")
    warning(msg_threshold)
  }
  if (!(V<5 & V>1)) {
    msg_folds <- paste("Number of folds:",V, "is either smaller than 1 or greater than 10")
    warning(msg_folds)
  }
  if (missing(SL.library)) {
    SL.library <- c("SL.mean", "SL.glm")
  }
  if (!is.data.frame(X)){
    X <- as.data.frame(X)
  }
  cvControl <- SuperLearner::SuperLearner.CV.control(V = V)
  method <- "method.NNLS"
  n <- nrow(X)
  objects <- vector("list", length(folds))
  for (ff in 1:length(folds)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(Xi[train_indices],
                                                data.frame(X, A = A)[train_indices, , drop = FALSE],
                                                SL.library = SL.library,
                                                family = binomial(),
                                                cvControl = cvControl, method = method)
  }
  nu <- function(a, x, ff) {
    if (length(a) != nrow(x)) stop("Length of 'a' must match number of rows in 'x'")
    out <- stats::predict(objects[[ff]],newdata = data.frame(cbind(x,A=a)))$pred
    out <- pmax(threshold, pmin(1 - threshold, out))
    return(out)
  }
  
  return(nu)
}

#' Estimate propensity score 
#'
#' This function trains the propensity score models 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data).
#' @param folds A list of cross-validation folds, typically created with \code{SuperLearner::CVFolds}. 
#' @param SL.library Vector of libraries for training SuperLearner.
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-3 by default). Constrains estimation to $[threshold, 1 - threshold]$.
#' @return A fold-specific function predicting propensity score given treatment (A) and covariates (X)
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' A <- stats::rbinom(100, 1, 0.5)
#' JFold <- 3
#' folds <- folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y, cvControl = SuperLearner::SuperLearner.CV.control(V = JFold, shuffle = TRUE))
#' SL.library <- c("SL.glm", "SL.mean")
#' prop_score <- estimate_ps(A, X, folds)
#' # Apply a function from the list to new data:
#' prop_score(1,X[1:5, ])  # Predict treated outcomes
#' }
#' @export
estimate_ps <- function(A, X, folds, SL.library, V = 2L, threshold = 1e-3) {
  if (!(threshold<5*1e-2 & threshold>0)) {
    msg_threshold <- paste("Threshold:", threshold, "is not small enough or not positive")
    warning(msg_threshold)
  }
  if (!(V<5 & V>1)) {
    msg_folds <- paste("Number of folds:",V, "is either smaller than 1 or greater than 10")
    warning(msg_folds)
  }
  if (missing(SL.library)) {
    SL.library <- c("SL.mean", "SL.glm")
  }
  if (!is.data.frame(X)){
    X <- as.data.frame(X)
  }
  #cvControl = list(V = 0)
  cvControl <- SuperLearner::SuperLearner.CV.control(V = V)
  method <- "method.CC_nloglik"
  n <- nrow(X)
  objects <- vector("list", length(folds))
  for (ff in 1:length(folds)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(A[train_indices],
                                                X[train_indices,],
                                                SL.library = SL.library,
                                                family = binomial(),
                                                method = method,
                                                cvControl = cvControl)
  }
  ps <- function(a, x, ff) {
    if (length(a) != nrow(x)) stop("Length of 'a' must match number of rows in 'x'")
    
    out <- rep(NA, length(a))
    if (sum(a == 1) > 0) {
      out[a == 1] <- stats::predict(objects[[ff]], newdata = x[a == 1, , drop = FALSE])$pred
    }
    if (sum(a == 0) > 0) {
      out[a == 0] <- 1-stats::predict(objects[[ff]], newdata = x[a == 0, , drop = FALSE])$pred
    }
    out <- pmax(threshold, pmin(1 - threshold, out))
    return(out)
  }
  
  return(ps)
}