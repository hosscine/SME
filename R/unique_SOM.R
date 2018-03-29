
#' R6 generator object of plane formed SOM.
#'
#' \code{planesom$new(dim, r, c, weights = NULL, alpha = 0.1, sigma = 1, collect.stats = F)}
#'
#' @export
#' @docType class
#' @format An \code{R6Class} generator object.
#'
panelsom <- R6Class(
  classname = "panelsom",
  inherit = tpsom,

  public = list(
    tprow = 4,
    tpcol = 8,

    initialize = function(dim, r = self$tprow, c = self$tpcol, weights = NULL,
                          neighbor = 1, alpha = 0.1, sigma = 1, collect.stats = F){

      super$initialize(dim = dim, skip.topology = T, weights=weights, neighbor = neighbor,
                       alpha = alpha, sigma = sigma, collect.stats = collect.stats)

      # instancing process
      if(!missing(r) || !missing(c)){
        self$setTopology((self$planeTopology(r,c)))
        self$tprow <- r
        self$tpcol <- c
      }

      # initializing process
      self$steps <- 0
      self$stats <- list()
      if(!is.null(self$weights.init)) self$weights <- self$weights.init
      else self$weights <- matrix(0,self$nnodes,dim)
    },

    planeTopology = function(row, col){
      tptmp <- matrix(1:(row * col), row, col)
      tp <- matrix(0, row + 2, col + 2)
      tp[2:(nrow(tp) - 1), 2:(ncol(tp) - 1)] <- tptmp
      ind <- expand.grid(2:(nrow(tp) - 1), 2:(ncol(tp) - 1))
      raw <- cbind(tp[rowMinus(ind, c(1, 0))], tp[rowMinus(ind, c(-1 , 0))],
                   tp[rowMinus(ind, c(0, 1))], tp[rowMinus(ind, c(0, -1))])
      lapply(1:nrow(raw), function(n){
        tmp <- raw[n,]
        tmp[tmp > 0]
      })
    }
  ))
