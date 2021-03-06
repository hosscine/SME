% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/som_generator.R
\name{setSom}
\alias{setSom}
\title{Generates various type of SOM with 20 nodes.}
\usage{
setSom(dim, topo.mode, weights = NULL, alpha = 0.1, sigma = 1,
  neighbor = 1, collect.stats = F)
}
\arguments{
\item{dim}{dimension of input vector.}

\item{topo.mode}{type of SOM topology. Choose from c(grid, plane, circle, sphere, torus, cylinder, cylinder.twocircle).}

\item{weights}{initial weights matrix.}

\item{alpha}{learning rate.}

\item{sigma}{variance of the weights propagation.}

\item{neighbor}{size of neighborhood kernel about the nodes. 1 indicates the neighbors of a node are its one adjacents.}

\item{collect.stats}{whether rocords learning history.}
}
\value{
SOM objsect.
}
\description{
Generates various type of SOM with 20 nodes.
}
\examples{
spsom <- setSom(3,"sphere")

}
