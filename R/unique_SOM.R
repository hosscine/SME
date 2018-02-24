#' create panel formed som
#' @export
#' @param dim dimension of input vector
#' @param r row of panel topology
#' @param c column of panel topology
setPanelSom <- function(dim,r,c) planesom$new(dim,r,c)

#' R6 generator object of plane formed SOM
#'
#' \code{planesom$new(dim, r, c, weights = NULL, alpha = 0.1, sigma = 1, collect.stats = F)}
#' @docType class
#' @export
#' @format An \code{R6Class} generator object.
planesom <- R6Class(
  classname = "planesom",
  inherit = "tpsom",

  public = list(
    tprow = 4,
    tpcol = 8,

    initialize = function(dim, r = self$tprow, c = self$tpcol, weights = NULL,
                          neighbor = 1, alpha = 0.1, sigma = 1, collect.stats = F){
      self$setTopology((planeTopology(r,c)))
      self$tprow <- r
      self$tpcol <- c

      self$dim <- dim
      if(!missing(weights)) self$weights <- weights
      # else self$weights <- matrix(runif(nnodes*dim),nnodes,dim)
      else self$weights <- matrix(0,nnodes,dim)

      self$neighbor.hop <- neighbor
      self$alpha <- 0.1
      self$sigma <- 1
      self$steps <- 0
      self$collect.stats <- collect.stats
    },

    planeTopology = function(r,c){
      tptmp <- matrix(1:(r*c),r,c)
      tp <- matrix(0,r+2,c+2)
      tp[2:(nrow(tp)-1),2:(ncol(tp)-1)] <- tptmp
      ind <- expand.grid(2:(nrow(tp)-1),2:(ncol(tp)-1))
      raw <- cbind(tp[rowMinus(ind,c(1,0))],tp[rowMinus(ind,c(-1,0))],
                   tp[rowMinus(ind,c(0,1))],tp[rowMinus(ind,c(0,-1))])
      lapply(1:nrow(raw),function(n){
        tmp <- raw[n,]
        tmp[tmp>0]
      })
    }
  ))

#' create a type of som
#' @export
#' @param dim Dimension of input vector.
#' @param topo.mode Id of typical som topology. Choose from c(grid, plane, circle, sphere, torus, cylinder, cylinder.twocircle).
#' @return Som objsect.
#' @examples spsom <- setSom(3,"sphere")
setSom <- function(dim, topo.mode){
  if(topo.mode == "grid")
    topo <- list(c(2,7,6),c(1,7,3),c(2,7,8,9,4),c(3,9,5),c(4,9,10),c(1,7,12,11),
                 c(1,2,3,8,12,6),c(3,9,14,13,12,7),c(3,4,5,10,14,8),c(5,9,14,15),
                 c(6,12,17,16),c(6,7,8,13,17,11),c(8,14,19,18,17,12),c(8,9,10,15,19,13),
                 c(10,20,19,14),c(11,17),c(11,12,13,18,16),c(17,13,19),c(18,13,14,15,20),
                 c(19,15))
  else if(topo.mode == "plane")
    topo <- list(c(2,6),c(1,3,7),c(2,4,8),c(3,5,9),c(4,10),c(1,7,11),c(2,6,8,12),c(3,7,9,13),
                 c(4,8,10,14),c(5,9,15),c(6,12,16),c(7,11,13,17),c(8,12,14,18),c(9,13,15,19),
                 c(10,14,20),c(11,17),c(12,16,18),c(13,17,19),c(14,18,20),c(15,19))
  else if(topo.mode == "circle")
    topo <- list(c(20,2),c(1,3),c(2,4),c(3,5),c(4,6),c(5,7),c(6,8),c(7,9),c(8,10),c(9,11),
                 c(10,12),c(11,13),c(12,14),c(13,15),c(14,16),c(15,17),c(16,18),c(17,19),
                 c(18,20),c(19,1))
  else if(topo.mode == "sphere")
    topo <- list(c(2,16,12),c(1,3,13),c(2,4,14),c(3,5,16),c(4,6,17),c(5,7,14),c(6,8,15),
                 c(7,9,17),c(8,10,18),c(9,11,15),c(10,12,13),c(1,11,18),c(2,11,19),
                 c(3,6,19),c(7,10,19),c(1,4,20),c(5,8,20),c(9,12,20),c(13,14,15),c(16,17,18))
  else if(topo.mode == "torus")
    topo <- list(c(2,4,13,16),c(1,3,5,14),c(2,6,15,16),c(1,5,7,17),c(2,4,6,8),c(3,5,9,17),
                 c(4,8,10,18),c(5,7,9,11),c(6,8,12,18),c(7,11,13,19),c(8,10,12,14),
                 c(9,11,16,19),c(1,10,14,20),c(2,11,13,15),c(3,12,14,20),c(1,3,17,20),
                 c(4,6,16,18),c(7,9,17,19),c(10,12,18,20),c(10,12,19,20))
  else if(topo.mode == "cylinder")
    topo <- list(c(2,5,6),c(1,3,7),c(2,4,8),c(3,5,9),c(1,4,10),
                 c(1,7,10,11),c(2,6,8,12),c(3,7,9,13),c(4,8,10,14),c(5,6,9,15),
                 c(6,12,15,16),c(7,11,13,17),c(8,12,14,18),c(9,13,15,19),c(10,11,14,20),
                 c(11,20,17),c(12,16,18),c(13,17,19),c(14,18,20),c(15,16,19))
  else if(topo.mode == "cylinder.twocircle")
    topo <- list(c(2,10,11),c(1,3,12),c(2,4,13),c(3,5,14),c(4,6,15),c(5,7,16),c(6,8,17),
                 c(7,9,18),c(8,10,19),c(9,1,20),c(1,12,20),c(2,11,13),c(3,12,14),
                 c(4,13,15),c(5,14,16),c(6,15,17),c(7,16,18),c(8,17,19),c(9,18,20),c(10,11,19))
  else stop("must chose topomode from c(grid, plane, circle, sphere, torus, cylinder, cylinder.twocircle)")

  som <- tpsom$new(dim = dim, topology = topo, nei = ifelse(topo.mode=="circle",2,1))
  return(som)
}

#' create graph of manifold from som
#' @export
#' @param som Corresponding som.
#' @param dim embedding dimension.
setGraph <- function(som,dim=NULL) tpgrp$new(som$adjacency,dim)
