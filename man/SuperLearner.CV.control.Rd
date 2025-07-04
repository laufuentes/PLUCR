% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross-fitting.R
\name{SuperLearner.CV.control}
\alias{SuperLearner.CV.control}
\title{SuperLearner.CV.control (from SuperLearner package)}
\usage{
SuperLearner.CV.control(
  V = 10L,
  stratifyCV = FALSE,
  shuffle = TRUE,
  validRows = NULL
)
}
\arguments{
\item{V}{Number of splits for the V-fold cross-validation step. The default is 10. In most cases, between 10 and 20 splits works well.}

\item{stratifyCV}{Should the data splits be stratified by a binary response? Attempts to maintain the same ratio in each training and validation sample.}

\item{shuffle}{Should the rows of X be shuffled before creating the splits.}

\item{validRows}{Use this to pass pre-specified rows for the sample splits. The length of the list should be V and each entry in the list should contain a vector with the row numbers of the corresponding validation sample.}
}
\value{
A list containing the control parameters.
}
\description{
SuperLearner.CV.control (from SuperLearner package)
}
