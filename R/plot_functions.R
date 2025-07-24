#' Plot Evolution of Objective Terms Across Lambda Values
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

#' Visualize Treatment Assignment Probability
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
    ggplot2::scale_color_continuous(limits = c(0, 1), oob = scales::squish) +
    ggplot2::labs(title = bquote(lambda == .(lambda) ~ "," ~ beta == .(beta))) +
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

#' Visualize the evolution of the RMSE and maximum RSE between consecutive iterations.  
#' 
#' Plots the evolution of the RMSE and RSE between iterative optimal solutions psi.  
#' 
#' @param intermediate_result File from Intermediate folder gathering results from alternated procedure. 
#' @param root.path Path to the folder where images are to be saved.
#' @param name A string to add to the end of filename. 
#' 
#' @return A message indicating that the image was saved.
#' @export
iterative_consecutive_metrics_plot <- function(intermediate_result, root.path, name){
  `%>%`<- magrittr::`%>%`
  psi_collection <- intermediate_result$psi_collection
  rmse_values <- numeric(ncol(psi_collection) - 1)
  max_rse_values <- numeric(ncol(psi_collection) - 1)
  
  for (i in seq_along(rmse_values)) {
    prev_psi <- psi_collection[,i]
    current_psi <- psi_collection[,i + 1]
    rmse_values[i] <- sqrt(mean((prev_psi - current_psi)^2))
    max_rse_values[i] <- max(sqrt((prev_psi - current_psi)^2))
  } 
  
  df <- tibble(
    iteration = 1:length(rmse_values),
    rmse_values = rmse_values, 
    max_rse_values = max_rse_values)
  
  rmse_plot <- ggplot2::ggplot(df, aes(x = iteration, y = rmse_values)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_point(color = "darkred") +
    ggplot2::labs(title = "RMSE Between Consecutive psi solutions",
      x = "Iteration",
      y = "RMSE") +
    ggplot2::theme_minimal()
  
  max_rse_plot <- ggplot2::ggplot(df, aes(x = iteration, y = max_rse_values)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_point(color = "darkred") +
    ggplot2::labs(title = "Max RSE Between Consecutive psi solutions",
         x = "Iteration",
         y = "Max RSE") +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(rmse_plot, filename = file.path(root.path, "Images",paste0("Iterative_RMSE_", name,".pdf")))
  ggplot2::ggsave(max_rse_plot, filename = file.path(root.path, "Images",paste0("Iterative_Max_RSE_", name,".pdf")))
  return("Images saved")
}

#' Visualize the evolution of the solutions psi obtained iteratively.  
#' 
#' Plot of the density of tje sigma_psi solutions obtained iteratively and the qq plot of the .  
#' 
#' @param intermediate_result File from Intermediate folder gathering results from alternated procedure. 
#' @param theta_opt The real numeric matrix (k x d). Each row is from oracular FW inner minimization, used to recover an extremal point for convex function construction.
#' @param theta_t An estimated numeric matrix (k x d). Each row is from T-learner FW inner minimization, used to recover an extremal point for convex function construction.
#' @param X_test A matrix of covariates from test data, of size n_test x d (input data).
#' @param beta A non-negative numeric scalar controlling the sharpness of the probability function (0.05 by default).
#' @param centered A logical value indicating whether to apply centering in \code{sigma_beta} (FALSE by default).
#' @param root.path Path to the folder where images are to be saved.
#' @param name A string to add to the end of filename. 
#' @return A message indicating that the image was saved.
#' @export
iterative_psi_evolution <- function(intermediate_result, theta_opt, theta_t, X_test, 
                                    beta=0.05, centered=FALSE, 
                                    root.path, name){
  `%>%`<- magrittr::`%>%`
  qq_plots <- NULL
  density_plots <- NULL
  ecdf_plots <- NULL
  
  # intermediate_result$theta_collection[[length(intermediate_result$theta_collection)+1]]<- intermediate_result$last_theta
  for(i in 1:length(intermediate_result$theta_collection)){
    theta_corr <- intermediate_result$theta_collection[[i]]
    current_psi <- make_psi(theta_corr)(X_test)
    
    df <- tibble::tibble(ora=PLUCR::sigma_beta(make_psi(theta_opt)(X_test), beta= 0.05, centered=centered),
                 tlearner=PLUCR::sigma_beta(make_psi(theta_t)(X_test), beta=0.05, centered=centered),
                 corr=PLUCR::sigma_beta(current_psi, beta=0.05, centered=centered))
    
    rmse_tlearner <- sqrt(mean((df$tlearner - df$ora)^2))
    rmse_corr <- sqrt(mean((df$corr - df$ora)^2))
    max_rse_tlearner <- max(sqrt((df$tlearner - df$ora)^2))
    max_rse_corr <- max(sqrt((df$corr - df$ora)^2))
    
    df_long_qq <- df %>%
      mutate(position=round(100 * (ecdf(ora)(ora)))) %>% 
      group_by(position) %>% slice_sample(n = 1)%>%
      pivot_longer(cols = c(tlearner, corr), names_to = "method", values_to = "value")
    
    p_qq <- ggplot(df_long_qq, aes(x = ora, y = value, color = method)) +
      geom_point(alpha = 0.6) + 
      geom_abline(intercept = 0, slope = 1, color = "red") +
      labs(
        subtitle = sprintf("RMSE: T-Learner = %.3f, Corrected = %.3f", rmse_tlearner, rmse_corr),
        x = "Oracle",
        y = "Estimated",
        color = "Method") +
      theme_minimal()
    
    qq_plots <- append(qq_plots, list(p_qq))
    
    df_long<- df %>%
      pivot_longer(cols = everything(), names_to = "method", values_to = "value")
    
    p_density <- ggplot(df_long, aes(x = value, color = method, fill = method)) +
      geom_density(alpha = 0.3) +
      labs(x = "sigma_beta Value",
           y = "Density") +
      theme_minimal()
    density_plots <- append(density_plots, list(p_density))
    
    p_ecdf <- ggplot(df_long, aes(x = value, color = method)) +
      stat_ecdf(geom = "step", direction = "hv") +  # 'hv' makes it ECQF-like
      coord_flip() +  # Flip axes to show quantile function (x = quantile)
      labs(
        x = "Quantile",
        y = "Value",
        color = "Method"
      ) +
      theme_minimal()
    ecdf_plots <- append(ecdf_plots, list(p_ecdf))
  }
  full_qq <- wrap_plots(qq_plots, ncol = 5) +
    plot_annotation(title = "QQ Plot: Oracle vs Estimated psi over Iterations")
  ggsave(filename = file.path(root.path, "Images", paste0("psi_evol_", name, ".pdf")), plot = full_qq, width = 19, height = 10)
  
  full_density <- wrap_plots(density_plots, ncol = 5) +
    plot_annotation(title = "Density Plot of sigma_beta Outputs")
  ggsave(filename = file.path(root.path, "Images", paste0("sigma_psi_density_", name, ".pdf")), plot = full_density, width = 19, height = 10)
  
  full_ecdf <- wrap_plots(ecdf_plots, ncol = 5) +
    plot_annotation(title = "Empirical Quantile Function (ECQF) of sigma_beta")
  ggsave(filename = file.path(root.path, "Images", paste0("sigma_psi_ecdf_", name, ".pdf")), plot = full_ecdf, width = 19, height = 10)
  
  return("Images saved")
}

#' synthetic_data_plot:: Plot Synthetic Data Setting
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
#' @return Saves a plot to "figures/synthetic_setting.pdf".
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
  df <- tidyr::expand_grid(x=seq(0,1,length.out=B), 
                    y=seq(0,1,length.out=B))
  df$delta_mu<-my_delta_Mu(df)
  df$delta_nu<-my_delta_Nu(df)
  df <- df %>% 
    tidyr::pivot_longer(cols = starts_with("delta"), 
                 names_to = 'What', 
                 values_to = 'Values')
  ggplot2::ggplot(df)+
    ggplot2::geom_raster(ggplot2::aes(x=x,y=y,fill=Values))+
    ggplot2::facet_grid(~What)+ 
    ggplot2::scale_fill_viridis_c(option = "magma")
  
  ggplot2::ggsave(file.path("root.path", "Images", paste0("Synthetic_data_plot_",name,".pdf")))
}



