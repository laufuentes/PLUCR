#' Plot Synthetic Data Setting
#'
#' Generates and saves a two-panel plot:
#' one showing the sign of the treatment effect (`delta_Mu`) and the other
#' visualizing the magnitude of selection effect (`delta_Nu`) across covariates X.1 and X.2.
#'
#' @param df_complete A data frame containing covariates prefixed with "X.".
#' @param delta_Mu A function that computes the treatment effect (mu difference) from covariates.
#' @param delta_Nu A function that computes the selection effect (nu difference) from covariates.
#'
#' @return Saves a plot to "figures/synthetic_setting.pdf".
#' @export
synthetic_setting_plot <- function(df_complete, delta_Mu, delta_Nu) {
  `%>%` <- magrittr::`%>%`
  df_complete$sign_delta_Mu <- as.factor(
    sign(
      delta_Mu(df_complete %>% dplyr::select(dplyr::starts_with("X.")))
    )
  )
  df_complete$delta_Nu <- delta_Nu(
    df_complete %>% dplyr::select(dplyr::starts_with("X."))
  )
  plot_Y_sign <- ggplot2::ggplot(df_complete, ggplot2::aes(x = df_complete$X.1, y = df_complete$X.2, color = df_complete$sign_delta_Mu)) +
    ggplot2::geom_point(alpha = 0.5)
  p_plot <- ggplot2::ggplot(df_complete, ggplot2::aes(x = df_complete$X.1, y = df_complete$X.2, color = (df_complete$delta_Nu))) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_gradient(low = "blue", high = "green")
  combined_plot <- gridExtra::grid.arrange(
    plot_Y_sign, p_plot,
    ncol = 1
  )
  ggplot2::ggsave(file.path("figures", "synthetic_setting.pdf"), combined_plot)
}

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
#'
#' @return A message indicating that the image was saved.
#' @export
visual_treatment_plot <- function(psi_X, lambda, beta, centered, Var_X_axis, Var_Y_axis, root.path) {
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
  
  ggsave(p, filename = file.path(root.path, "Images","Treatment_assignment"))
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


#' Visualize the evolution of the correction term for mu and nu. 
#' 
#' Plots the correction terms for mu and nu over iterations of the alternated procedure 
#' 
#' @param intermediate_result File from Intermediate folder gathering results from alternated procedure. 
#' @param root.path Path to the folder where images are to be saved.
#' @param name A string to add to the end of filename. 
#' 
#' @return A message indicating that the image was saved.
#' @export
iterative_bias_plot <- function(intermediate_result, root.path, name){
  
  bias_evol <- data.frame(
    iterations = seq(1, length(intermediate_result$correction_term_mu), by = 1),
    mu_correction_term = t(intermediate_result$correction_term_mu),
    nu_correction_term = t(intermediate_result$correction_term_nu)
  )
  
  bias_evol_long <- bias_evol %>%
    pivot_longer(
      cols = c("mu_correction_term", "nu_correction_term"), 
      names_to = "bias_type", 
      values_to = "bias")
  
  bias_plot <- ggplot(bias_evol_long, aes(x = iterations, y = bias, color=bias_type)) +
    geom_line() +
    geom_point() +
    theme_minimal()
  
  ggsave(bias_plot, filename = file.path(root.path, "Images",paste0("Iterative_bias_",name,".pdf")))
  
  term1_df <- as.data.frame(intermediate_result$term1)  # samples x iterations
  term2_df <- as.data.frame(intermediate_result$term2)
  
  term1_stats <- data.frame(
    iterations = 1:ncol(term1_df),
    mean = colMeans(term1_df),
    variance = apply(term1_df, 2, var),
    term = "term1"
  )
  
  term2_stats <- data.frame(
    iterations = 1:ncol(term2_df),
    mean = colMeans(term2_df),
    variance = apply(term2_df, 2, var),
    term = "term2"
  )
  
  mean_plot_term1 <- ggplot(term1_stats, aes(x = iterations, y = mean)) +
    geom_line(color = "blue") + geom_point(color = "blue") +
    labs(title = "Mean of term1 over iterations") +
    theme_minimal()
  
  # term1 variance plot
  var_plot_term1 <- ggplot(term1_stats, aes(x = iterations, y = variance)) +
    geom_line(color = "blue") + geom_point(color = "blue") +
    labs(title = "Variance of term1 over iterations") +
    theme_minimal()
  
  # term2 mean plot
  mean_plot_term2 <- ggplot(term2_stats, aes(x = iterations, y = mean)) +
    geom_line(color = "red") + geom_point(color = "red") +
    labs(title = "Mean of term2 over iterations") +
    theme_minimal()
  
  # term2 variance plot
  var_plot_term2 <- ggplot(term2_stats, aes(x = iterations, y = variance)) +
    geom_line(color = "red") + geom_point(color = "red") +
    labs(title = "Variance of term2 over iterations") +
    theme_minimal()
  
  wrap_plots(list(mean_plot_term1, mean_plot_term2,var_plot_term1,var_plot_term2), ncol = 2) 
  ggsave(filename = file.path(root.path, "Images", paste0("Bias_terms_", name, ".pdf")), width = 19, height = 10)
  
  return("Images saved")
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
  
  rmse_plot <- ggplot(df, aes(x = iteration, y = rmse_values)) +
    geom_line(color = "steelblue") +
    geom_point(color = "darkred") +
    labs(title = "RMSE Between Consecutive psi solutions",
      x = "Iteration",
      y = "RMSE") +
    theme_minimal()
  
  max_rse_plot <- ggplot(df, aes(x = iteration, y = max_rse_values)) +
    geom_line(color = "steelblue") +
    geom_point(color = "darkred") +
    labs(title = "Max RSE Between Consecutive psi solutions",
         x = "Iteration",
         y = "Max RSE") +
    theme_minimal()
  
  ggsave(rmse_plot, filename = file.path(root.path, "Images",paste0("Iterative_RMSE_", name,".pdf")))
  ggsave(max_rse_plot, filename = file.path(root.path, "Images",paste0("Iterative_Max_RSE_", name,".pdf")))
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
  qq_plots <- NULL
  density_plots <- NULL
  ecdf_plots <- NULL
  
  # intermediate_result$theta_collection[[length(intermediate_result$theta_collection)+1]]<- intermediate_result$last_theta
  for(i in 1:length(intermediate_result$theta_collection)){
    theta_corr <- intermediate_result$theta_collection[[i]]
    current_psi <- make_psi(theta_corr)(X_test)
    
    df <- tibble(ora=PLUCR::sigma_beta(make_psi(theta_opt)(X_test), beta= 0.05, centered=centered),
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
        title = "Comparison of ORA with T-Learner and Corrected",
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
      labs(title = "Density Plot of sigma_beta Outputs",
           x = "sigma_beta Value",
           y = "Density") +
      theme_minimal()
    density_plots <- append(density_plots, list(p_density))
    
    p_ecdf <- ggplot(df_long, aes(x = value, color = method)) +
      stat_ecdf(geom = "step", direction = "hv") +  # 'hv' makes it ECQF-like
      coord_flip() +  # Flip axes to show quantile function (x = quantile)
      labs(
        title = "Empirical Quantile Function (ECQF)",
        x = "Quantile",
        y = "Value",
        color = "Method"
      ) +
      theme_minimal()
    ecdf_plots <- append(ecdf_plots, list(p_ecdf))
  }
  wrap_plots(qq_plots, ncol = 5)
  ggsave(filename = file.path(root.path,"Images",paste0("psi_evol_",name,".pdf")), width = 19, height = 10)
  
  wrap_plots(density_plots, ncol = 5)
  ggsave(filename = file.path(root.path,"Images",paste0("sigma_psi_density_",name,".pdf")), width = 19, height = 10)
  
  wrap_plots(ecdf_plots, ncol = 5)
  ggsave(filename = file.path(root.path,"Images",paste0("sigma_psi_ecdf_",name,".pdf")), width = 19, height = 10)
  return("Images saved")
}

#' Plot Optimal Policy Across Discrete Lambda Values
#'
#' Generates multiple plots of treatment assignment probability for a range of lambda values
#' and combines them with a shared legend. Also saves plots for the initial and optimal lambdas.
#'
#' param results A list or data frame containing `optimal_x`, `lambda`, `beta`, etc.
#' param idx_opt The index of the optimal lambda in the results.
#' param df A data frame with covariates, including X.1 and X.2.
#' param type_simu A string used to name the output folder (e.g., "oracular", "estimated").
#' param centered Logical; if TRUE, uses centered sigma transformation.
#'
#' return Saves plots in "figures/<type_simu>/".
#' export
#geom_points_fct <- function(results, idx_opt, df, type_simu, centered){
#  lambda_discr <- which(results$beta==results$beta[idx_opt])[as.integer(seq(1, length(which(results$beta==results$beta[idx_opt])), length.out=10))]
#  
#  plots <- lapply(lambda_discr, function(x) gamma_plot_funct(results$optimal_x[[x]], results$lambda[[x]], results$beta[[x]], df, centered))
#  plots_no_legend <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
  
#  legend <- cowplot::get_legend(plots[[1]])
#  combined_plots <- cowplot::plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
#  final_plot <- cowplot::plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
  # Display the final plot
#  print(final_plot)
#  ggplot2::ggsave(file.path("figures",type_simu,"optimal_solution_multiple_lambdas.pdf"),final_plot, width = 10, height = 6)
  
#  plot_none<-gamma_plot_funct(results$optimal_x[[1]], results$lambda[[1]], results$beta[[1]], df, centered)+ggplot2::theme(legend.position = "none")
#  plot_max <- gamma_plot_funct(results$optimal_x[[idx_opt]], results$lambda[[idx_opt]], results$beta[[idx_opt]], df, centered)+ggplot2::theme(legend.position = "none")

#  opt_plots <- cowplot::plot_grid(plot_none,plot_max, ncol=2)
#  ggplot2::ggsave(file.path("figures",type_simu,"optimal_solution_optimal_lambda.pdf"),opt_plots, width = 8, height = 4)
#}