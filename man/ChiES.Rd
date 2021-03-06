% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChiES.R
\name{ChiES}
\alias{ChiES}
\title{Compute power for Chi Square Based on Effect Size
Takes phi, degrees of freedom, and a range of sameple sizes. Alpha is .05 by default, alterative values may be entered by user}
\usage{
ChiES(phi, df, nlow, nhigh, by = 1, alpha = 0.05)
}
\arguments{
\item{phi}{phi coefficient (effect size for 2x2)}

\item{df}{degrees of freedom}

\item{nlow}{starting sample size}

\item{nhigh}{ending sample size}

\item{by}{Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)}

\item{alpha}{Type I error (default is .05)}
}
\value{
Power for Chi Square Based on Effect Size
}
\description{
Compute power for Chi Square Based on Effect Size
Takes phi, degrees of freedom, and a range of sameple sizes. Alpha is .05 by default, alterative values may be entered by user
}
