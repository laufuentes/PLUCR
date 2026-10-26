#' Visualize treatment assignment probability
#'
#' Plots the smoothed treatment assignment probability over covariates Var_X_axis and Var_Y_axis
#'
#' @param psi_X A numeric vector with values in the range `[-1,1]` representing the output of \code{psi} at fixed X.
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