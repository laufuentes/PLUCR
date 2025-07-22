#' SL.grf
#'
#' This function trains conditional mean of a target variable for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param Y outcome variable
#' @param X training dataframe
#' @param newX test dataframe
#' @param family gaussian or binomial 
#' @param obsWeights observation-level weights 
#' @param ... not used 
#' @return a list containing the predictions and the fitted object
#' 
#' @export
SL.grf <- function(Y, X, newX, family, obsWeights, ...) {
  if (is.matrix(X)) X <- as.data.frame(X)
  if (is.matrix(newX)) newX <- as.data.frame(newX)
  
  is_binary <- setequal(unique(Y), c(0, 1))
  
  fit.grf <- if (is_binary) {
    fit.grf <- grf::probability_forest(Y = as.factor(Y), X = X)
  } else {
    fit.grf <- grf::regression_forest(Y = Y, X = X)
  }
  
  pred_output <- predict(fit.grf, newdata = newX)
  
  # Extract predictions based on the type of forest
  pred <- if (is_binary) {
    # If it's a probability_forest, predictions is a matrix with P(Y=0) and P(Y=1)
    # SuperLearner (binomial family) expects P(Y=1), which is usually the second column
    # or explicitly named "1".
    if (is.matrix(pred_output$predictions) && ncol(pred_output$predictions) == 2) {
      as.numeric(pred_output$predictions[, "1"]) # Prefer named column for robustness
    } else {
      # Fallback in case it's not a 2-column matrix (e.g., older GRF version)
      as.numeric(pred_output$predictions)
    }
  } else {
    # If it's a regression_forest, predictions is already a vector.
    # No need to subscript with [,1]
    as.numeric(pred_output$predictions)
  }
  
  # The crucial check that helps identify length mismatches
  if (length(pred) != nrow(newX)) {
    stop("Prediction length mismatch in SL.grf: Expected ", nrow(newX), " predictions, got ", length(pred))
  }
  
  fit <- list(object = fit.grf)
  class(fit) <- "SL.grf"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' predict.SL.grf
#'
#' This function trains conditional mean of primary outcome models for treated and control groups 
#' using `SuperLearner`, applying cross-validation to compute out-of-fold estimates.
#'
#' @param object SL.grf object
#' @param newdata dataframe to generate predictions
#' @param ... not used 
#' @return the requested predictions 
#' 
#' @export
predict.SL.grf <- function (object, newdata, ...) {
  if (is.matrix(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  # Get predictions from the underlying grf forest
  pred_output <- predict(object = object$object, newdata = newdata)
  
  # Determine if the original model was a probability forest or regression forest
  # This can be inferred from the class of the stored grf object.
  if (inherits(object$object, "probability_forest")) {
    # If it was a probability_forest, extract P(Y=1)
    if (is.matrix(pred_output$predictions) && ncol(pred_output$predictions) == 2) {
      pred <- as.numeric(pred_output$predictions[, "1"]) # Prefer named column
    } else {
      # Fallback
      pred <- as.numeric(pred_output$predictions)
    }
  } else { # Assume it's a regression_forest or other continuous outcome model
    # If it was a regression_forest, predictions is already a vector
    pred <- as.numeric(pred_output$predictions)
  }
  return(pred)
}
