
#' R6 generator object of plane formed SOM.
#'
#' \code{planesom$new(dim, r, c, weights = NULL, alpha = 0.1, sigma = 1, collect.stats = F)}
#'
#' @importFrom myfs padding.matrix
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
        self$setTopology((private$planeTopology(r,c)))
        self$tprow <- r
        self$tpcol <- c
      }

      # initializing process
      self$steps <- 0
      self$stats <- list()
      if(!is.null(self$weights.init)) self$weights <- self$weights.init
      else self$weights <- matrix(0,self$nnodes,dim)
    }),

  private = list(
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

    })
)


#' R6 generator object of torus formed SOM.
#'
#' \code{torussom$new(dim, r, c, weights = NULL, alpha = 0.1, sigma = 1, collect.stats = F)}
#'
#' @export
#' @docType class
#' @format An \code{R6Class} generator object.
#'
torussom <- R6Class(
  classname = "torussom",
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
        self$setTopology((private$torusAdjacent(r,c)))
        self$tprow <- r
        self$tpcol <- c
      }

      # initializing process
      self$steps <- 0
      self$stats <- list()
      if(!is.null(self$weights.init)) self$weights <- self$weights.init
      else self$weights <- matrix(0,self$nnodes,dim)
    }),

  private = list(
    torusAdjacent = function(nr, nc){
      plane <- matrix(1:(nr * nc), nr, nc)
      plane.torus <- padding.matrix(plane, size = 1, replace = 0)
      adjacent <- matrix(0, nr * nc, nr * nc)

      # fill borders of padding matrix by opposit side of the borders
      plane.torus[1, 2:(nc + 1)] <- plane[nr,]
      plane.torus[2:(nr + 1), 1] <- plane[, nc]
      plane.torus[nr + 2, 2:(nc + 1)] <- plane[1,]
      plane.torus[2:(nr + 1), nc + 2] <- plane[,1]

      # find adjacents from torusian plane
      adjacent[cbind(1:(nr * nc), plane.torus[1:nr, 2:(nc + 1)] %>% as.numeric)] <- 1
      adjacent[cbind(1:(nr * nc), plane.torus[2:(nr + 1), 1:nc] %>% as.numeric)] <- 1
      adjacent[cbind(1:(nr * nc), plane.torus[3:(nr + 2), 2:(nc + 1)] %>% as.numeric)] <- 1
      adjacent[cbind(1:(nr * nc), plane.torus[2:(nr + 1), 3:(nc + 2)] %>% as.numeric)] <- 1

      return(adjacent)
    })
)
