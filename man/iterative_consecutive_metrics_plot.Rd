% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{iterative_consecutive_metrics_plot}
\alias{iterative_consecutive_metrics_plot}
\title{Visualize the evolution of the RMSE and maximum RSE between consecutive iterations.}
\usage{
iterative_consecutive_metrics_plot(intermediate_result, root.path, name)
}
\arguments{
\item{intermediate_result}{File from Intermediate folder gathering results from alternated procedure.}

\item{root.path}{Path to the folder where images are to be saved.}

\item{name}{A string to add to the end of filename.}
}
\value{
A message indicating that the image was saved.
}
\description{
Plots the evolution of the RMSE and RSE between iterative optimal solutions psi.
}
