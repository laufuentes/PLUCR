#' SL.grf
#'
#' This function trains conditional mean of primary outcome models for treated and control groups 
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
SL.grf <- function (Y, X, newX, family, obsWeights, ...) {
    if (is.matrix(X)) {
        X <- as.data.frame(X)
    }
    fit.grf <- grf::regression_forest(Y = Y, X = X)
    if (is.matrix(newX)) {
        newX <- as.data.frame(newX)
    }
    pred <- as.vector(predict(fit.grf, newdata = newX)$predictions)
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
    pred <- as.vector(predict(object = object$object, newdata = newdata)$predictions)
    pred
}
