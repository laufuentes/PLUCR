#' Plot evolution of objective terms across lambda values
#'
#' Visualizes how risk, constraint, and overall objective evolve with respect to different lambda values.
#' Includes smooth loess trends and confidence intervals.
#'
#' @param results A data frame or tibble containing `lambda`, `risk`, `constraint`, `obj`, and `beta` columns.
#' @param type_simu A string indicating the simulation type (e.g., "oracular" or "empirical").
#' @param beta_opt The optimal non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#'
#' @return Saves a plot to "figures/<type_simu>/lambda_evol.pdf".
#' @export
lambda_evol <- function(results, type_simu, beta_opt) {
  results_lambda <- results[which(results$beta == beta_opt), ]
  res <- as.data.frame(results_lambda)
  # lambda plot
  lambda_evol_plot <- ggplot2::ggplot(res, ggplot2::aes(x = res$lambda)) +
    ggplot2::geom_point(ggplot2::aes(y = res$risk, color = "Risk")) +
    ggplot2::geom_point(ggplot2::aes(y = res$constraint, color = "Constraint")) +
    #geom_point(aes(y = policy_value, color = "Policy value")) +
    ggplot2::geom_point(ggplot2::aes(y = res$obj, color = "L")) +
    # Risk: with colored SE and no line
    ggplot2::geom_smooth(ggplot2::aes(y = res$risk, color = "Risk", fill = "Risk"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    # Constraint: with colored SE and no line
    ggplot2::geom_smooth(ggplot2::aes(y = res$constraint, color = "Constraint", fill = "Constraint"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    ggplot2::geom_smooth(ggplot2::aes(y = res$obj, color = "L", fill = "L"),
                         method = "loess", se = TRUE, span = 0.75, size = 0, alpha = 0.2
    ) +
    ggplot2::labs(
      title = expression("Evolution of optimal solution " * psi[lambda] * " with respect to " * lambda),
      subtitle = if (type_simu == "oracular") {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[0]](psi, lambda) * ": " * psi %in% Psi), "}"))
      } else {
        expression(psi[lambda] == argmin~bgroup("{", list(L[P[n]](psi, lambda) * ": " * psi %in% Psi), "}"))
      },
      x = expression(lambda),
      y = "Values"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(fill = "none")
  
  ggplot2::ggsave(file.path("figures", type_simu, "lambda_evol.pdf"), lambda_evol_plot)
}

#' Visualize treatment assignment probability
#'
#' Plots the smoothed treatment assignment probability over covariates Var_X_axis and Var_Y_axis
#'
#' @param psi_X A numeric vector with values in the range \code{[-1, 1]} representing the output of \code{psi} at fixed X.
#' @param lambda Regularization scalar for the constraint (display purposes).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param Var_X_axis Vector of covariates values for x-axis (length n).
#' @param Var_Y_axis Vector of covariates values for y-axis (length n).
#' @param root.path Path to the folder where images are to be saved.
#' @param name A string to add to the end of filename. 
#'
#' @return A message indicating that the image was saved.
#' @export
visual_treatment_plot <- function(psi_X, lambda, beta, centered, Var_X_axis, Var_Y_axis, root.path, name) {
  policy <- sigma_beta(psi_X, beta, centered)
  
  df <- data.frame(
    x = Var_X_axis,
    y = Var_Y_axis,
    treat_proba = policy
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = treat_proba)) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::scale_color_viridis_c(option = "magma", limits = c(0, 1), oob = scales::squish) +
    ggplot2::labs(
      title = bquote(lambda == .(lambda) ~ "," ~ beta == .(beta)),
      color = "Treatment\nProbability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  
  ggsave(p, filename = file.path(root.path, "Images",paste0("Treatment_assignment_",name,".pdf")))
  return("Image saved")
} 

#' Visualize the true value versus the estimated version. 
#' 
#' Plots the density of the true values for treated and control and the estimated versions. 
#' 
#' @param X A matrix of covariates of size n x d (input data).
#' @param Var1 Vector of treated variable values (length n).
#' @param Var0 Vector of control variable values  (length n).
#' @param Var_learner A function predicting a variable (Y, Xi or propensity score) given treatment (A) and covariates (X).
#' @param root.path Path to the folder where images are to be saved.
#' @param variable_name A string precising the name of the variable ("Y", "Xi", or "prop_score")
#' 
#' @return A message indicating that the image was saved.
#' @export
plot_nuisance <- function(X, Var1, Var0, Var_learner, root.path, variable_name){
  `%>%`<- magrittr::`%>%`
  pred1 <- Var_learner(rep(1,nrow(X)),X)
  pred0 <- Var_learner(rep(0,nrow(X)),X)
  
  df <- tibble(Oracular.1 = Var1, 
               Estimated.1 = pred1, 
               Oracular.0 = Var0, 
               Estimated.0 = pred0)
  
  df_long <- df %>% 
    tidyr::pivot_longer(
      cols = everything(),
      names_to = c("type", "outcome"),
      names_pattern = "(Oracular|Estimated)\\.([01])",
      values_to = "value"
    )
  
  n <- ggplot(df_long, aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ outcome, scales = "free") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(title = "Density Comparison: Oracular vs Estimated Potential Outcomes",
         x = "Outcome Value", y = "Density", fill = "Type") +
    ggplot2::theme_minimal()
  
  ggsave(n, filename = file.path(root.path, "Images",paste0("Nuisance_", variable_name, ".pdf")))
  return("Image saved")
}

#' Plot synthetic data setting
#'
#' Generates and saves a two-panel plot:
#' one showing the sign of the treatment effect (`delta_Mu`) and the other
#' visualizing the magnitude of selection effect (`delta_Nu`) across covariates X.1 and X.2.
#'
#' @param delta_Mu A function that computes the treatment effect (mu difference) from covariates.
#' @param delta_Nu A function that computes the selection effect (nu difference) from covariates.
#' @param B Integer, number of Monte Carlo repetitions (1e4 by default).
#' @param root.path Path to the folder where images are to be saved.
#' @param name A string to add to the end of filename. 
#' @return Saves a plot to "Images/synthetic_setting.pdf".
#' @export
synthetic_data_plot <-function(delta_Mu, delta_Nu, B=1e2, root.path, name){
  `%>%`<- magrittr::`%>%`
  # Add proper xlab and ylab for representation
  vars_mu <- attr(delta_Mu, "vars")
  vars_nu <- attr(delta_Nu, "vars")
  my_delta_Mu <- function(df){
    df <- as.matrix(df)
    mat <- matrix(0, nrow=nrow(df), ncol = max(vars_mu))
    for(i in 1:2){
      mat[,vars_mu[i]]<- df[,i]
    }
    return(delta_Mu(mat))
  } 
  my_delta_Nu <- function(df){
    df <- as.matrix(df)
    mat <- matrix(0, nrow=nrow(df), ncol = max(vars_nu))
    for(i in 1:2){
      mat[,vars_nu[i]]<- df[,i]
    }
    return(delta_Nu(mat))
  } 
  
  x_col <- paste0("X", vars_nu[1])
  y_col <- paste0("X", vars_nu[2])
  df <- tidyr::expand_grid(x=seq(0,1,length.out=B), 
                    y=seq(0,1,length.out=B))
  df$delta_mu<-my_delta_Mu(df)
  df$delta_nu<-my_delta_Nu(df)
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = starts_with("delta"),
      names_to = "What",
      values_to = "Values"
    )
  
  df_long$What <- factor(df_long$What, levels = c("delta_mu", "delta_nu"))
  
  # Set the labels as strings that label_parsed can parse
  levels(df_long$What) <- c("Delta * mu[0](X)", "Delta * nu[0](X)")
  
  p<- ggplot(df_long) +
    geom_raster(aes(x = x, y = y, fill = Values)) +
    facet_grid(~What, labeller = label_parsed) +  # label_parsed will parse the factor levels
    scale_fill_viridis_c(option = "magma", limits = c(-1, 1)) +
    labs(x = "X1", y = "X2") +
    theme_minimal()
  
  ggplot2::ggsave(file.path(root.path, "Images", paste0("Synthetic_data_plot_",name,".pdf")), width = 5, height = 4)
}

#' Plot metric values for comparison
#'
#' Creates a comparison plot of different metrics across different treatment rule estimation methods.
#'
#' The function takes a data frame with method names and corresponding policy values, constraints, etc. (typically including `Theta_0`, `Theta_naive`, `Theta_final`, and `Theta_oracular`) and generates a visual comparison (e.g., bar plot or point plot).
#'
#' #' @param data A tibble or data frame with two columns:
#' \describe{
#'   \item{method}{A character vector indicating the estimation method (e.g., "Theta_0", "Theta_naive", etc.).}
#'   \item{metric1}{A numeric vector containing the corresponding metric1 values for each method.}
#'   \item{metric2}{A numeric vector containing the corresponding metric2 values for each method.}}
#' @param metrics A vector containing the metrics to be represented. 
#' @param techniques A vector containing the names of the techniques to be represented. 
#' @param root.path Path to the folder where images are to be saved.
#' @return Saves a plot to "Images/"Comparison_techniques.pdf".
#' @export
#' 
plot_metric_comparison <- function(data, 
                                   metrics = c("policy_value", "constraint"),
                                   techniques = NULL,
                                   root.path) {
  # Pivot to long format: one row per (method, metric)
  data_long <- data %>%
    pivot_longer(cols = all_of(metrics), 
                 names_to = "metric", 
                 values_to = "value")
  
  # Basic grouped bar plot
  plot <- ggplot2::ggplot(data_long, aes(x = metric, y = value, fill = method)) +
    ggplot2::geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    ggplot2::labs(
      title = "Comparison of Methods",
      x = "Metric",
      y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.title = element_blank())
  
  # Save the plot
  filename <- paste0("Comparison_techniques.pdf")
  ggplot2::ggsave(plot, filename = file.path(root.path, "Images", filename), width = 6, height = 4)
  
  return(plot)
}