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
  
  ggplot2::ggsave(p, filename = file.path(root.path, "Images",paste0("Treatment_assignment_",name,".pdf")))
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
synthetic_data_plot <- function(delta_Mu, delta_Nu, B = 1e2, root.path, name) {
  `%>%` <- magrittr::`%>%`
  
  vars_mu <- attr(delta_Mu, "vars")
  vars_nu <- attr(delta_Nu, "vars")
  
  my_delta_Mu <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = max(vars_mu))
    for (i in 1:2) mat[, vars_mu[i]] <- df[, i]
    delta_Mu(mat)
  }
  my_delta_Nu <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = max(vars_nu))
    for (i in 1:2) mat[, vars_nu[i]] <- df[, i]
    delta_Nu(mat)
  }
  
  df <- tidyr::expand_grid(
    x = seq(0, 1, length.out = B),
    y = seq(0, 1, length.out = B)
  )
  df$delta_mu <- my_delta_Mu(df)
  df$delta_nu <- my_delta_Nu(df)
  
  p_mu <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = delta_mu)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(option = "magma", limits = c(-1, 1)) +
    ggplot2::labs(
      title = expression(Delta * mu[0](X)),
      x = paste0("X", vars_mu[1]),
      y = paste0("X", vars_mu[2]),
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")   # <- no legend here
  
  p_nu <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = delta_nu)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(option = "magma", limits = c(-1, 1)) +
    ggplot2::labs(
      title = expression(Delta * nu[0](X)),
      x = paste0("X", vars_nu[1]),
      y = paste0("X", vars_nu[2]),
      fill = NULL
    ) +
    ggplot2::theme_minimal()                   # <- legend stays here
  
  p <- gridExtra::grid.arrange(p_mu, p_nu, ncol = 2)
  
  ggplot2::ggsave(
    filename = file.path(root.path, "Images", paste0("Synthetic_data_plot_", name, ".pdf")),
    plot = p,
    width = 8, height = 4
  )
}


#' Plot realistic data setting
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
plot_realistic <- function(delta_Mu, delta_Nu, B = 100, root.path, name) {
  `%>%` <- magrittr::`%>%`
  vars_mu <- attr(delta_Mu, "vars")
  vars_nu <- attr(delta_Nu, "vars")
  exp<- generate_realistic_data(2)
   
  
  my_delta_Mu <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = max(vars_mu))
    for (i in seq_along(vars_mu)) {
      mat[, vars_mu[i]] <- df[, i]
    }
    attr(mat, "max_Y")<- unique(exp[[1]]$max_Y)
    attr(mat, "min_Y")<-unique(exp[[1]]$min_Y)
    delta_Mu(mat)
  }
  
  my_delta_Nu <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = max(vars_nu))
    for (i in seq_along(vars_nu)) {
      mat[, vars_nu[i]] <- df[, i]
    }
    attr(mat, "max_Y")<- unique(exp[[1]]$max_Y)
    attr(mat, "min_Y")<-unique(exp[[1]]$min_Y)
    delta_Nu(mat)
  }
  
  # Generate grid for mu (continuous) and nu (binary) separately
  # vars_mu: continuous ranges
  x_mu <- as.integer(seq(16, 65, length.out = B))  # for column 1
  y_mu <- seq(-5, 5, length.out = B)   # for column 4 (adjust as needed)
  
  # vars_nu: binary ranges
  x_nu <- c(0,1)  # column 2
  y_nu <- c(0,1)  # column 3
  
  df_mu <- tidyr::expand_grid(x = x_mu, y = y_mu)
  df_mu$delta_mu <- my_delta_Mu(df_mu)
  
  df_nu <- tidyr::expand_grid(x = x_nu, y = y_nu)[-2,]
  df_nu$delta_nu <- my_delta_Nu(df_nu)
  
  # Now column names are the same: x, y, delta_*
  df_mu_long <- df_mu %>% tidyr::pivot_longer(cols = delta_mu,
                                              names_to = "What", values_to = "Values")
  df_nu_long <- df_nu %>% tidyr::pivot_longer(cols = delta_nu,
                                              names_to = "What", values_to = "Values")
  
  df_long <- rbind(df_mu_long, df_nu_long)
  df_long$What <- factor(df_long$What, levels = c("delta_mu", "delta_nu"))
  levels(df_long$What) <- c("Delta * mu[0](X)", "Delta * nu[0](X)")
  
  # Plot
  p_mu <- ggplot2::ggplot(df_mu_long) +
    ggplot2::geom_raster(ggplot2::aes(x = x,y = y,fill = Values))+
    ggplot2::scale_fill_viridis_c(option = "magma", limits=c(-1,1)) +
    ggplot2::labs(x = "X.1", y = "X.4") +
    ggplot2::theme(legend.position = "none")
  
  p_nu <- ggplot2::ggplot(df_nu_long) +
    ggplot2::geom_raster(ggplot2::aes(x = x,y = y,fill = Values))+
    ggplot2::scale_fill_viridis_c(option = "magma", limits=c(-1,1)) +
    ggplot2::labs(x = "X.2", y = "X.3") +
    ggplot2::theme_minimal()
  
  p <- gridExtra::grid.arrange(p_mu, p_nu, ncol = 2)
  
  ggplot2::ggsave(file.path(root.path, "Images", paste0("Synthetic_data_plot_", name, ".pdf")), plot = p)
}



#' Plot metric values for comparison
#'
#' Creates a comparison plot of different metrics across different treatment rule estimation methods.
#'
#' The function takes a data frame with method names and corresponding policy values, constraints, etc. (typically including `Theta_0`, `Theta_naive`, `Theta_final`, and `Theta_oracular`) and generates a visual comparison (e.g., bar plot or point plot).
#'
#' @param data A tibble or data frame with two columns:
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
  `%>%`<- magrittr::`%>%`
  data_long <- data %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(metrics), 
                 names_to = "metric", 
                 values_to = "value")
  
  # Basic grouped bar plot
  plot <- ggplot2::ggplot(data_long, ggplot2::aes(x = metric, y = value, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    ggplot2::labs(
      title = "Comparison of Methods",
      x = "Metric",
      y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
  # Save the plot
  filename <- paste0("Comparison_techniques.pdf")
  ggplot2::ggsave(plot, filename = file.path(root.path, "Images", filename), width = 6, height = 4)
  
  return(plot)
}