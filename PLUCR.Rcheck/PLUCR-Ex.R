pkgname <- "PLUCR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PLUCR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("HX")
### * HX

flush(stderr()); flush(stdout())

### Name: HX
### Title: Compute the Inverse Propensity Score Weight (IPW)
### Aliases: HX

### ** Examples

# Example usage:
prop_model <- function(A, X) { 0.5 }  # Constant propensity score for illustration
HX(1, data.frame(x1 = 1, x2 = 2), prop_model)




cleanEx()
nameEx("Optimization_Estimation")
### * Optimization_Estimation

flush(stderr()); flush(stdout())

### Name: Optimization_Estimation
### Title: Iterative optimization procedure
### Aliases: Optimization_Estimation

### ** Examples

# (Requires user-defined functions: mu0, nu0, prop_score, FW, make_psi, sigma_beta, update_mu_XA, update_nu_XA)
# Optimization_Estimation(mu0, nu0, prop_score, df, 
                          lambda=1, alpha=0.1, precision=0.025, beta=0.05, 
                         centered=TRUE, folder="path/to/folder", prefix="run1")




cleanEx()
nameEx("delta_mu_constant")
### * delta_mu_constant

flush(stderr()); flush(stdout())

### Name: delta_mu_constant
### Title: Constant Conditional Average Treatment Effect estimator for Y
### Aliases: delta_mu_constant

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_constant(X)



cleanEx()
nameEx("delta_mu_linear")
### * delta_mu_linear

flush(stderr()); flush(stdout())

### Name: delta_mu_linear
### Title: Linear-shaped Conditional Average Treatment Effect estimator for
###   Y
### Aliases: delta_mu_linear

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_linear(X)



cleanEx()
nameEx("delta_mu_mix")
### * delta_mu_mix

flush(stderr()); flush(stdout())

### Name: delta_mu_mix
### Title: Mixed-shape Conditional Average Treatment Effect estimator for Y
### Aliases: delta_mu_mix

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_mix(X)



cleanEx()
nameEx("delta_mu_null")
### * delta_mu_null

flush(stderr()); flush(stdout())

### Name: delta_mu_null
### Title: Null Conditional Average Treatment Effect estimator for Y
### Aliases: delta_mu_null

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_null(X)



cleanEx()
nameEx("delta_mu_realistic")
### * delta_mu_realistic

flush(stderr()); flush(stdout())

### Name: delta_mu_realistic
### Title: Realistic Conditional Average Treatment Effect estimator for Y
### Aliases: delta_mu_realistic

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_realistic(X)



cleanEx()
nameEx("delta_mu_threshold")
### * delta_mu_threshold

flush(stderr()); flush(stdout())

### Name: delta_mu_threshold
### Title: Thresholded-shaped Conditional Average Treatment Effect
###   estimator for Y
### Aliases: delta_mu_threshold

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_mu_threshold(X)



cleanEx()
nameEx("delta_nu_linear")
### * delta_nu_linear

flush(stderr()); flush(stdout())

### Name: delta_nu_linear
### Title: Linear-shaped Conditional Average Treatment Effect estimator for
###   Xi
### Aliases: delta_nu_linear

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_nu_linear(X)



cleanEx()
nameEx("delta_nu_mix")
### * delta_nu_mix

flush(stderr()); flush(stdout())

### Name: delta_nu_mix
### Title: Mixed-shaped Conditional Average Treatment Effect estimator for
###   Xi
### Aliases: delta_nu_mix

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_nu_mix(X)



cleanEx()
nameEx("delta_nu_realistic")
### * delta_nu_realistic

flush(stderr()); flush(stdout())

### Name: delta_nu_realistic
### Title: Realistic Conditional Average Treatment Effect estimator for Xi
### Aliases: delta_nu_realistic

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_nu(X)



cleanEx()
nameEx("delta_nu_satisfied")
### * delta_nu_satisfied

flush(stderr()); flush(stdout())

### Name: delta_nu_satisfied
### Title: Computes the difference in expected outcomes under treatment and
###   control.
### Aliases: delta_nu_satisfied

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_nu_satisfied(X)



cleanEx()
nameEx("delta_nu_threshold")
### * delta_nu_threshold

flush(stderr()); flush(stdout())

### Name: delta_nu_threshold
### Title: Thresholded Conditional Average Treatment Effect estimator for
###   Xi
### Aliases: delta_nu_threshold

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
delta_nu_threshold(X)



cleanEx()
nameEx("estimate_mu")
### * estimate_mu

flush(stderr()); flush(stdout())

### Name: estimate_mu
### Title: Estimate mu
### Aliases: estimate_mu

### ** Examples

## Not run: 
##D set.seed(123)
##D X <- matrix(stats::rnorm(100 * 5), ncol = 5)
##D A <- stats::rbinom(100, 1, 0.5)
##D Y <- stats::runif(100)
##D JFold <- 3
##D folds <- SuperLearner::CVFolds(
##D   n, 
##D   id = NULL,
##D   Y = Y, 
##D   cvControl = SuperLearner::SuperLearner.CV.control(
##D     V = JFold, 
##D     shuffle = TRUE
##D   )
##D )
##D SL.library <- c("SL.glm", "SL.mean", "SL.ranger")
##D mu_functions <- estimate_mu(Y, A, X, folds)
##D # Apply a function from the list to new data:
##D mu_functions(1,X[1:5, ],1)  # Predict treated outcomes
## End(Not run)



cleanEx()
nameEx("estimate_nu")
### * estimate_nu

flush(stderr()); flush(stdout())

### Name: estimate_nu
### Title: Estimate nu
### Aliases: estimate_nu

### ** Examples

## Not run: 
##D set.seed(123)
##D X <- matrix(stats::rnorm(100 * 5), ncol = 5)
##D A <- stats::rbinom(100, 1, 0.5)
##D Xi <- stats::rbinom(100, n=1, p=0.5)
##D JFold <- 3
##D folds <- SuperLearner::CVFolds(
##D   n, 
##D   id = NULL,
##D   Y = Y, 
##D   cvControl = SuperLearner::SuperLearner.CV.control(
##D     V = JFold, 
##D     shuffle = TRUE
##D   )
##D )
##D SL.library <- c("SL.glm", "SL.mean", "SL.ranger")
##D nu_functions <- estimate_nu(Xi, A, X, folds)
##D # Apply a function from the list to new data:
##D nu_functions(1,X[1:5, ],1)  # Predict treated outcomes
## End(Not run)



cleanEx()
nameEx("estimate_ps")
### * estimate_ps

flush(stderr()); flush(stdout())

### Name: estimate_ps
### Title: Estimate propensity score
### Aliases: estimate_ps

### ** Examples

## Not run: 
##D set.seed(123)
##D X <- matrix(stats::rnorm(100 * 5), ncol = 5)
##D A <- stats::rbinom(100, 1, 0.5)
##D JFold <- 3
##D folds <- SuperLearner::CVFolds(
##D   n, 
##D   id = NULL,
##D   Y = Y, 
##D   cvControl = SuperLearner::SuperLearner.CV.control(
##D     V = JFold, 
##D     shuffle = TRUE
##D   )
##D )
##D SL.library <- c("SL.glm", "SL.mean")
##D prop_score <- estimate_ps(A, X, folds)
##D # Apply a function from the list to new data:
##D prop_score(1,X[1:5, ])  # Predict treated outcomes
## End(Not run)



cleanEx()
nameEx("generate_data")
### * generate_data

flush(stderr()); flush(stdout())

### Name: generate_data
### Title: Synthetic data generator and functions generator
### Aliases: generate_data

### ** Examples

data <- data_gen(100)
head(data[[1]])  # complete data
head(data[[2]])  # observed data



cleanEx()
nameEx("generate_realistic_data")
### * generate_realistic_data

flush(stderr()); flush(stdout())

### Name: generate_realistic_data
### Title: Realistic synthetic data generator and functions generator
### Aliases: generate_realistic_data

### ** Examples

data <- data_gen(100)
head(data[[1]])  # complete data
head(data[[2]])  # observed data



cleanEx()
nameEx("model_Xi_linear")
### * model_Xi_linear

flush(stderr()); flush(stdout())

### Name: model_Xi_linear
### Title: Linear treatment effect on Xi Component Function
### Aliases: model_Xi_linear

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
model_Xi_linear(X)



cleanEx()
nameEx("model_Xi_mix")
### * model_Xi_mix

flush(stderr()); flush(stdout())

### Name: model_Xi_mix
### Title: Mixed treatment effect on Xi component function
### Aliases: model_Xi_mix

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
model_Xi_mix(X)



cleanEx()
nameEx("model_Xi_realistic")
### * model_Xi_realistic

flush(stderr()); flush(stdout())

### Name: model_Xi_realistic
### Title: Realistic treatment effect on Xi Component Function
### Aliases: model_Xi_realistic

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
model_Xi_linear(X)



cleanEx()
nameEx("model_Xi_satisfied")
### * model_Xi_satisfied

flush(stderr()); flush(stdout())

### Name: model_Xi_satisfied
### Title: Low treatment effect on Xi
### Aliases: model_Xi_satisfied

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
model_Xi_satisfied(X)



cleanEx()
nameEx("model_Xi_threshold")
### * model_Xi_threshold

flush(stderr()); flush(stdout())

### Name: model_Xi_threshold
### Title: Thresholded treatment effect on Xi component function
### Aliases: model_Xi_threshold

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
model_Xi_threshold(X)



cleanEx()
nameEx("model_Y_constant")
### * model_Y_constant

flush(stderr()); flush(stdout())

### Name: model_Y_constant
### Title: Constant treatment effect on Y component function
### Aliases: model_Y_constant

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_constant_TE(X, A)



cleanEx()
nameEx("model_Y_linear")
### * model_Y_linear

flush(stderr()); flush(stdout())

### Name: model_Y_linear
### Title: Linear treatment effect on Y component function
### Aliases: model_Y_linear

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_linear(X, A)



cleanEx()
nameEx("model_Y_mix")
### * model_Y_mix

flush(stderr()); flush(stdout())

### Name: model_Y_mix
### Title: Mixed treatment effect on Y component function
### Aliases: model_Y_mix

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_mix(X, A)



cleanEx()
nameEx("model_Y_null")
### * model_Y_null

flush(stderr()); flush(stdout())

### Name: model_Y_null
### Title: No treatment effect on Y component function
### Aliases: model_Y_null

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_null_TE(X, A)



cleanEx()
nameEx("model_Y_realistic")
### * model_Y_realistic

flush(stderr()); flush(stdout())

### Name: model_Y_realistic
### Title: Realistic treatment effect on Y component function
### Aliases: model_Y_realistic

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_linear(X, A)



cleanEx()
nameEx("model_Y_threshold")
### * model_Y_threshold

flush(stderr()); flush(stdout())

### Name: model_Y_threshold
### Title: Thresholded treatment effect on Y component function
### Aliases: model_Y_threshold

### ** Examples

X <- matrix(stats::runif(10*2), 10, 2)
A <- rep(1, 10)
model_Y_threshold(X, A)



cleanEx()
nameEx("phi")
### * phi

flush(stderr()); flush(stdout())

### Name: phi
### Title: Normalize a Matrix by Column Min-Max Scaling
### Aliases: phi

### ** Examples

# Example matrix
X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
X_norm <- phi(X)
X_norm
attributes(X_norm)$min_X
attributes(X_norm)$max_X




cleanEx()
nameEx("phi_inv")
### * phi_inv

flush(stderr()); flush(stdout())

### Name: phi_inv
### Title: Inverse Min-Max Normalization
### Aliases: phi_inv

### ** Examples

# Normalize and then invert
X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
X_norm <- phi(X)
X_recovered <- phi_inv(X_norm)
all.equal(X, X_recovered)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
