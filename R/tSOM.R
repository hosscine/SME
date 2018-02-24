#' topological SOM
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An R6Class generator object.
tpsom <-
  R6Class(
    classname = "tpsom", inherit = tpgrp,
    public = list(
      # Public Fields ----------------------------------------------------------
      # nnodes = "numeric",
      # topology = "list",
      # adjacency = "matrix",
      # hop = "matrix",
      #
      # dim = "numeric",
      # weights = "matrix",

      neighbor.hop = 1,
      alpha = 0.1,
      sigma = 1,
      steps = 0,
      stats = list(),
      collect.stats = F,

      # Public Methods ----------------------------------------------------------
      initialize = function(dim, topology, weights = NULL,
                            neighbor = 1, alpha = 0.1, sigma = 1, collect.stats = F){
        self$setTopology(topology)

        self$dim <- dim
        if(!missing(weights)) self$weights <- weights
        # else self$weights <- matrix(runif(nnodes*dim),nnodes,dim)
        else self$weights <- matrix(0,self$nnodes,dim)

        self$neighbor.hop <- neighbor
        self$alpha <- 0.1
        self$sigma <- 1
        self$steps <- 0
        self$collect.stats <- collect.stats
      },

      updateWeights = function(x, win)
        for (i in self$calcNeighbor(win))
          self$weights[i,] <- self$weights[i,] + self$hci(win,i) * (x - self$weights[i,]),

      hci = function(c, i) self$alpha*exp(-self$calcHop(c,i)/2/self$sigma^2),

      collectStats = function(x, win){
        if(!self$collect.stats)
          self$stats <- append(self$stats,list(c(step=self$steps,
                                                 input=x,
                                                 win=win,
                                                 error=vnorm(x-self$weights[win,]),
                                                 distortion=self$calcDistortion(x,win)
          )))
      },

      calcDistortion = function(x, win)
        sum(sapply(self$calcNeighbor(win),
                   function(n) self$hci(win,n) * vnorm(x-self$weights[n,]))),

      step = function(x){
        self$steps <- steps + 1
        winner <- self$calcWinner(x)
        self$collectStats(x,winner)
        self$updateWeights(x,winner)
        return(winner)
      },

      batchStep = function(X,t=1000)
        hoge <- apply(X[round(runif(t)*nrow(X)),],1,function(x)step(x))
    )
  )
