
#' Sets learning parameters for SOM
#'
#' @param som target som.
#' @param alpha learning rate.
#' @param sigma variance of the weights propagation.
#' @param neighbor size of neighborhood kernel about the nodes. 1 indicates the neighbors of a node are its one adjacents.
#'
#' @export
#'
setLearningPar <- function(som, alpha = 0.1, sigma = 1, neighbor = 2){
  som$alpha <- alpha
  som$sigma <- sigma
  som$neighbor.hop <- neighbor
}
