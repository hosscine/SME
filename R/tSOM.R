#' R6 generator object of topological SOM.
#'
#' \code{tpsom$new(dim, topology, weights = NULL, neighbor = 1, alpha = 0.1, sigma = 1, collect.stats = F)}
#'
#' @importFrom stats prcomp
#' @importFrom rgl points3d
#' @importFrom rgl segments3d
#'
#' @export
#' @docType class
#' @format An \code{R6Class} generator object.
#'
tpsom <-
  R6Class(
    classname = "tpsom", inherit = tpgrp,
    public = list(
      # Public Fields ----------------------------------------------------------
      neighbor.hop = 1,
      alpha = 0.1,
      sigma = 1,
      collect.stats = F,

      steps = 0,
      stats = list(),

      weights.init = NULL,

      # Public Methods ----------------------------------------------------------
      initialize = function(dim, topology, weights=NULL,
                            neighbor = 1, alpha = 0.1, sigma = 1, collect.stats = F){
        # check args
        if(missing(topology) && is.null(self$adjacency))
          stop("to instancing the class, you must specify a topology.")
        if(missing(dim) && is.null(self$dim))
          stop("to instancing the class, you must specify a dimension.")

        # instancing process
        if(!missing(topology)) self$setTopology(topology)
        if(!missing(dim)) self$dim <- dim
        if(!missing(neighbor)) self$neighbor.hop <- neighbor
        if(!missing(alpha)) self$alpha <- 0.1
        if(!missing(sigma)) self$sigma <- 1
        if(!missing(collect.stats)) self$collect.stats <- collect.stats
        if(!missing(weights)) self$weights.init <- weights

        self$steps <- 0
        self$stats <- list()
        if(!is.null(self$weights.init)) self$weights <- self$weights.init
        else self$weights <- matrix(0,self$nnodes,dim)

      },

      updateWeights = function(x, win)
        for (i in self$calcNeighbor(win))
          self$weights[i,] <- self$weights[i,] + self$hci(win,i) * (x - self$weights[i,]),

      hci = function(c, i) self$alpha*exp(-self$calcHop(c,i)/2/self$sigma^2),

      collectStats = function(x, win){
        if(self$collect.stats)
          self$stats <- append(self$stats,list(c(step=self$steps,
                                                 input=x,
                                                 win=win,
                                                 error=vnorm(x-self$weights[win,]),
                                                 distortion=self$calcDistortion(x,win)
          )))
      },

      plot = function(X, ..., as.grp = F,
                      dim = min(sapply(1:self$nnodes,function(n)length(self$calcNeighbor(n))))-1){
        if(missing(X) || as.grp){
          super$plot(..., as.grp = T, dim = dim)
          cat("plot by graph mode.\n
              if you want to plot by data mode, add argment X and set as.grp = F.")
        }
        else {
          if(self$dim == 2){
            elp <- overwriteEllipsis(..., xlab = "", ylab = "", x = X, col = 1)
            do.call(plot,elp)
            points(self$weights, ...)
          }
          else if(self$dim == 3){
            elp <- overwriteEllipsis(..., xlab = "", ylab = "", zlab = "",
                                     x = X, col = 1, size = 5)
            do.call(plot3d,elp)
            points3d(self$weights, ...)
          }
          else{
            pca <- prcomp(rbind(X, self$weights))[[5]]

            # plot X
            elp <- overwriteEllipsis(..., xlab = "", ylab = "", zlab = "",
                                     x = pca[1:nrow(X),], size = 5,
                                     col = rainbow(self$nnodes)[apply(X, 1, self$calcWinner)])
            do.call(plot3d, elp)

            # plot SOM
            rgl::points3d(pca[(nrow(X)+1):nrow(pca),], size = 3)
            for(i in 1:self$nnodes){
              nei <- self$calcNeighbor(i, neighbor.hop = 1)
              rgl::segments3d(x = pca[replace(rep(i, 2 * length(nei)), 2 * 1:length(nei), nei) + nrow(X),])
            }
            rgl::text3d(pca[(nrow(X) + 1):nrow(pca),], texts = 1:self$nnodes, adj = 1)
          }
        }
      },

      calcDistortion = function(x, win)
        sum(sapply(self$calcNeighbor(win),
                   function(n) self$hci(win,n) * vnorm(x-self$weights[n,]))),

      step = function(x){
        self$steps <- self$steps + 1
        winner <- self$calcWinner(x)
        self$collectStats(x,winner)
        self$updateWeights(x,winner)
        return(winner)
      },

      batchStep = function(X,t=1000)
        hoge <- apply(X[round(runif(t)*nrow(X)),],1,function(x)self$step(x))
    )
  )
