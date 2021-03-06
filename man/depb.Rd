% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depb.R
\name{depb}
\alias{depb}
\title{Power for Comparing Dependent Coefficients in Multiple Regression with Two or Three Predictors
Requires correlatiosn between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)}
\usage{
depb(ry1, ry2, ry3 = NULL, r12, r13 = NULL, r23 = NULL, n = NULL,
  alpha = 0.05)
}
\arguments{
\item{ry1}{Correlation between DV (y) and first predictor (1)}

\item{ry2}{Correlation between DV (y) and second predictor (2)}

\item{ry3}{Correlation between DV (y) and third predictor (3)}

\item{r12}{Correlation between first (1) and second predictor (2)}

\item{r13}{Correlation between first (1) and third predictor (3)}

\item{r23}{Correlation between second (2) and third predictor (3)}

\item{n}{Sample size}

\item{alpha}{Type I error (default is .05)}

\item{rep}{number of replications (default is 10000)}
}
\value{
Power for Comparing Dependent Coefficients in Multiple Regression with Two or Three Predictors
}
\description{
Power for Comparing Dependent Coefficients in Multiple Regression with Two or Three Predictors
Requires correlatiosn between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
}
