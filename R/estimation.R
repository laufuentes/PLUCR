#' Estimate mu 
#'
#' This function trains conditional mean of primary outcome models for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in `[0,1]`).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data in `[0,1]`).
#' @param folds A list of cross-validation folds (e.g., a list of indices for each fold). 
#' @param SL.library Vector of libraries for training SuperLearner (c("SL.glm", "SL.mean") by default).
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-2 by default). Constrains estimation to `[threshold, 1 - threshold]`.
#' @return A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X)
#' @export
estimate_mu <- function(Y, A, X, folds, SL.library = c("SL.glm", "SL.mean"), V = 2L, threshold = 1e-2) {
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
  for (ff in c(1,3)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(Y[train_indices],
                                                data.frame(X, A = A)[train_indices, , drop = FALSE],
                                                SL.library = SL.library,
                                                family = stats::gaussian(),
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
#' @param X A matrix or data frame of covariates of size n x d (input data in `[0,1]`).
#' @param folds A list of cross-validation folds (e.g., a list of indices for each fold).
#' @param SL.library Vector of libraries for training SuperLearner (c("SL.glm", "SL.mean") by default).
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-2 by default). Constrains estimation to `[threshold, 1 - threshold]`.
#' @return A fold-specific function predicting adverse event outcome (Xi) given treatment (A) and covariates (X)
#' @export
estimate_nu <- function(Xi, A, X, folds, SL.library=c("SL.glm", "SL.mean"), V = 2L, threshold = 1e-2) {
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
  for (ff in c(1,3)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(Xi[train_indices],
                                                data.frame(X, A = A)[train_indices, , drop = FALSE],
                                                SL.library = SL.library,
                                                family = stats::binomial(),
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
#' @param X A matrix or data frame of covariates of size n x d (input data in `[0,1]`).
#' @param folds A list of cross-validation folds (e.g., a list of indices for each fold).
#' @param SL.library Vector of libraries for training SuperLearner (c("SL.glm", "SL.mean") by default).
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @param threshold A numeric scalar that sets the minimum allowed value for upper and lower bound estimations (1e-2 by default). Constrains estimation to `[threshold, 1 - threshold]`.
#' @return A fold-specific function predicting propensity score given treatment (A) and covariates (X)
#' @export
estimate_ps <- function(A, X, folds, SL.library=c("SL.glm", "SL.mean"), V = 2L, threshold = 1e-2) {
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
  method <- "method.CC_nloglik"
  n <- nrow(X)
  objects <- vector("list", length(folds))
  for (ff in c(1,3)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(A[train_indices],
                                                X[train_indices,],
                                                SL.library = SL.library,
                                                family = stats::binomial(),
                                                method = method,
                                                cvControl = cvControl)
  }
  ps <- function(a, x, ff) {
    if (length(a) != nrow(x)) stop("Length of 'a' must match number of rows in 'x'")
    
    out <- rep(NA, length(a))
    preds <- stats::predict(objects[[ff]], newdata = x[, drop = FALSE])$pred
    if (sum(a == 1) > 0) {
      out[a == 1] <- preds[a==1]
    }
    if (sum(a == 0) > 0) {
      out[a == 0] <- 1- preds[a==0]
    }
    out <- pmax(threshold, pmin(1 - threshold, out))
    return(out)
  }
  
  return(ps)
}


#' Estimate real-valued mu 
#'
#' This function trains conditional mean of primary outcome models for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param Y A numeric vector or matrix of length n representing primary outcomes (in R).
#' @param A A binary vector or matrix of length n indicating treatment assignment (0 or 1).
#' @param X A matrix or data frame of covariates of size n x d (input data in R).
#' @param folds A list of cross-validation folds (e.g., a list of indices for each fold).
#' @param SL.library Vector of libraries for training SuperLearner (c("SL.glm", "SL.mean") by default).
#' @param V Number of folds inside the SuperLearner (2L by default).
#' @return A fold-specific function predicting primary outcome (Y) given treatment (A) and covariates (X)
#' @export
estimate_real_valued_mu <- function(Y, A, X, folds, SL.library=c("SL.glm", "SL.mean"), V = 2L) {
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
  for (ff in c(1,3)) {
    train_indices <- folds[[ff]]
    objects[[ff]] <- SuperLearner::SuperLearner(Y[train_indices],
                                                data.frame(X, A = A)[train_indices, , drop = FALSE],
                                                SL.library = SL.library,
                                                family = stats::gaussian(),
                                                cvControl = cvControl, method = method)
  }
  mu <- function(a, x, ff) {
    if (length(a) != nrow(x)) stop("Length of 'a' must match number of rows in 'x'")
    out <- rep(NA, length(a))
    out <- stats::predict(objects[[ff]],newdata = data.frame(x,A=a))$pred
    return(out)
  }
  return(mu)
}