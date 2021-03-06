% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/som_generator.R
\name{setPlaneSom}
\alias{setPlaneSom}
\title{Creates plane formed SOM.}
\usage{
setPlaneSom(dim, r, c, ...)
}
\arguments{
\item{dim}{dimension of input vector.}

\item{r}{row of plane topology.}

\item{c}{column of plane topology.}

\item{...}{SOM generater's parameters such as weights = NULL, alpha = 0.1, sigma = 1, neighbor = 1, collect.stats = F.}
}
\value{
plane SOM object.
}
\description{
Creates plane formed SOM.
}
