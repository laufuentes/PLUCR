#' Generate Psi Function
#'
#' Constructs a convex combination function \code{psi} based on the sequence of solutions
#' obtained from the Frank-Wolfe (FW) algorithm. Each new solution \code{theta} contributes  
#' to \code{psi} in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
#'
#' @param Theta A numeric matrix (k x d). Each row \code{theta} is from FW inner minimization, used to recover an extremal point for convex function construction.
#'
#' @return A function \code{psi} that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]},
#'         using a convex combination of past \code{theta} solutions.
#' @export
make_psi <- function(Theta) {
  # Theta: k x d
  # gamma: real
  if(FALSE){
    ## ------------
    ## lazy version
    ## ------------
    psi_iterative <- function(x) {
      # x: n x d
      K <- nrow(Theta)
      Theta_x <- x %*% t(Theta)
      psi_x <- 2*expit(Theta_x[,1])-1
      weights <- rep(0,K)
      weights[1] <- 1
      # Theta_x: n x K
      for (k in 1:(ncol(Theta_x)-1)){
        gamma <- 2/(2+k)
        weights[1:k] <- weights[1:k]*(1-gamma)
        weights[k+1] <- gamma
        psi_x <- (1-gamma) * psi_x + gamma * (2 * expit(Theta_x[, k+1])-1)
      }
      return(psi_x)
    }
  }
  ## ----------------
  ## Smarter version
  ## ----------------
  psi <- function(x){
    k <- nrow(Theta)
    Gamma <- matrix((2 / ((k + 1) * k)) * (1:k), ncol = 1)
    sigma_theta_x <- 2 * expit(x %*% t(Theta)) - 1 # only use the first k thetas.
    psi_x <- sigma_theta_x %*% Gamma
    return(psi_x)
  }
  return(psi)
}

#' Stochastic Gradient Descent (SGD)
#'
#' Performs stochastic gradient descent to optimize the parameters.
#'
#' @param theta_current A numeric matrix of size 1 x d (intialization for parameter to estimate).
#' @param psi A function that takes an input \code{X} and returns a numeric vector with values in the range \code{[-1, 1]}.
#' @param X A matrix of covariates of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the contrast between primary outcomes.
#' @param delta_Nu A function of \code{X} that determines the contrast between adverse event outcomes.
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param batch_prop Proportion of data in a batch (by default 1/3).
#' @param max_iter Maximum number of iterations in the SGD (by default 1e3).
#' @param tol Tolerance parameter (by default 1e-3).
#' @param lr Learning rate parameter (by default 1e-2).
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{FALSE}.
#'
#' @return A numeric matrix of size 1 x d (optimized parameters).
#' @export
SGD <- function(theta_current, psi, X, delta_Mu, delta_Nu, lambda, alpha=0.1, beta=0.05, centered=FALSE,
                batch_prop=1/3, max_iter=1e3, tol=1e-3, lr=1e-2, verbose=FALSE){
  if(!(1/5 <= batch_prop & 4/5 >= batch_prop)){
    warning("Argument batch_prop in call to SGD is either small or large (that is, not in [1/5,4/5]).\n")
  }
  if(max_iter <= 10){
    warning("Argument max_iter in call to SGD is small (smaller than 10).\n")
  }
  if(tol <= 1e-5 | tol >= 1e-1){
    warning("Argument tol in call to SGD is either very small or very large (that is, not in [1e-5,1e-1]).\n")
  }
  if(lr <= 1e-4 | lr >= 1e-1){
    warning("Argument lr in call to SGD is either very small or very large (that is, not in [1e-4,1e-1]).\n")
  }
  
  n <- nrow(X)
  batch_size <- n*batch_prop
  
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  LprimeX <-  grad_Lagrangian_p(psi, X, delta_Mu, delta_Nu, lambda, alpha, beta, centered)
  for(i in 1:max_iter){
    s <- sample.int(n, batch_size)
    x <- X[s,]

    theta_x <- x %*% t(theta_current)
    expit_theta_x <- expit(theta_x)
    expit_diff <- 2 * expit_theta_x * (1 - expit_theta_x)
    
    Lprime <-LprimeX[s]
    dL_dtheta <- t(t(x) %*% (expit_diff * Lprime))
    theta_current <- theta_current - lr * dL_dtheta
    
    if (verbose && i %% 500 == 0) {
            theta_X <- X %*% t(theta_current)
            expit_theta_X_full <- expit(theta_X)
            expit_Diff <- 2 * expit_theta_X_full * (1 - expit_theta_X_full)

            whole_grad <- t(t(X) %*% (expit_Diff * LprimeX))

            if (mean(whole_grad) < tol) {
                break
            }
            value <- mean(LprimeX * (2 * expit_theta_X_full - 1))
            msg <- sprintf("\tSGD: iteration %i, value %f", i, value)
            message(msg)}
}
    return(theta_current)
}

#' @rdname SGD
#' 
#' @export
SGD_X <- function(theta_current, psi_X, X, delta_Mu_X, delta_Nu_X, lambda, alpha=0.1, beta=0.05, centered=FALSE,
                batch_prop=1/3, max_iter=1e3, tol=1e-3, lr=1e-2, verbose){
  if(!(1/5 <= batch_prop & 4/5 >= batch_prop)){
    warning("Argument batch_prop in call to SGD is either small or large (that is, not in [1/5,4/5]).\n")
  }
  if(max_iter <= 10){
    warning("Argument max_iter in call to SGD is small (smaller than 10).\n")
  }
  if(tol <= 1e-5 | tol >= 1e-1){
    warning("Argument tol in call to SGD is either very small or very large (that is, not in [1e-5,1e-1]).\n")
  }
  if(lr <= 1e-4 | lr >= 1e-1){
    warning("Argument lr in call to SGD is either very small or very large (that is, not in [1e-4,1e-1]).\n")
  }
  
  n <- nrow(X)
  batch_size <- n*batch_prop
  
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  LprimeX <-  grad_Lagrangian_p(psi_X, delta_Mu_X, delta_Nu_X, lambda, alpha, beta, centered)
  for(i in 1:max_iter){
    s <- sample.int(n, batch_size)
    x <- X[s,]
    
    theta_x <- x %*% t(theta_current)
    expit_theta_x <- expit(theta_x)
    expit_diff <- 2 * expit_theta_x * (1 - expit_theta_x)
    
    Lprime <- LprimeX[s]
    dL_dtheta <- t(t(x) %*% (expit_diff * Lprime))
    theta_current <- theta_current - lr * dL_dtheta
    
    if (verbose && i %% 500 == 0) {
      theta_X <- X %*% t(theta_current)
      expit_theta_X_full <- expit(theta_X)
      expit_Diff <- 2 * expit_theta_X_full * (1 - expit_theta_X_full)
      
      whole_grad <- t(t(X) %*% (expit_Diff * LprimeX))
      
      if (mean(whole_grad) < tol) {
        break
      }
      value <- mean(LprimeX * (2 * expit_theta_X_full - 1))
      msg <- sprintf("SGD: iteration %i, value %f", i, value)
      message(msg)}
  }
  return(theta_current)
}

#' Frank-Wolfe Algorithm
#'
#' Implements the Frank-Wolfe optimization algorithm to iteratively refine a convex  
#' combination function \code{psi}. At each iteration, a new solution \code{theta}  
#' is computed via stochastic gradient descent (SGD) and added to the convex combination  
#' in the form \eqn{2 \cdot \text{expit}(X \theta) - 1}.
#'
#' @param X A matrix of covariates of size n x d (input data).
#' @param delta_Mu A function of \code{X} that determines the contrast between primary outcomes.
#' @param delta_Nu A function of \code{X} that determines the contrast between adverse event outcomes.
#' @param lambda A non-negative numeric scalar controlling the penalty for violating the constraint.
#' @param alpha A numeric scalar representing the constraint tolerance (in [0,1/2], 0.1 by default).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical (FALSE by default) indicating whether to center the policy.
#' @param precision A numeric scalar defining the desired convergence precision (0.05 by default). The number of Frank-Wolfe iterations (K) is inversely proportional to this value, calculated as 1/precision.
#' @param verbose A logical value indicating whether to print progress updates. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the optimized parameter \code{theta},  
#'         where each row represents the k-th \code{theta} solution at iteration \code{k}.
#' @export
FW <- function(X, delta_Mu, delta_Nu, lambda, alpha=0.1, beta=0.05, centered=FALSE, precision=0.05, verbose=TRUE) {
    K <- as.integer(1/precision)
    d <- ncol(X)
    theta_init <- matrix(stats::runif(d, -5, 5), ncol=d, nrow=1)
    theta <- theta_init

    for (k in 0:K){
      if (k==1){theta <- matrix(theta[2,], nrow=1, ncol=d)}

      psi <- make_psi(theta)
        
        if (verbose && k %% 20 == 0) {
            msg <- sprintf("FW: iteration %i, value %f", k, Lagrangian_p(psi, X, delta_Mu, delta_Nu, lambda=lambda, alpha=alpha, beta=beta, centered=centered))
            message(msg)
        }
        theta_opt <- SGD(theta_current=theta_init, psi, X, delta_Mu, delta_Nu,
                         lambda, alpha, beta, centered, verbose=(verbose && k %% 10 == 0))
        
        theta <- rbind(theta, theta_opt)
    }
    return(theta)
}

#' @rdname FW
#' 
#' @export
FW_X <- function(X, delta_Mu, delta_Nu, lambda, alpha=0.1, beta=0.05, centered=FALSE, precision=0.05, verbose=TRUE) {
  K <- as.integer(1/precision)
  tol <- 1e-5
  d <- ncol(X)
  theta <- matrix(stats::runif(d, -5, 5), ncol=d, nrow=1)

  for (k in 0:K){
    if (k==1){theta <- matrix(theta[2,], nrow=1, ncol=d)}
    
    psi <- make_psi(theta)
    psi_X<- psi(X)
    delta_Mu_X <- delta_Mu(X)
    delta_Nu_X <- delta_Nu(X)
    if (verbose && k %% 20 == 0) {
      msg <- sprintf("FW: iteration %i, value %f", k, Lagrangian_p(psi, X, delta_Mu, delta_Nu, lambda, alpha, beta, centered))
      message(msg)
    }
    theta_opt <- SGD_X(theta_current=theta, psi_X, X, delta_Mu_X, delta_Nu_X,
                     lambda, alpha, beta, centered, verbose=(verbose && k %% 10 == 0))
    
    theta <- rbind(theta, theta_opt)
  }
  return(theta)
}