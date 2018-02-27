#' R6 generator object of topological graph.
#'
#' \code{tpgrp$new(topology, dim = NULL, weights = NULL)}
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom pipeR %>>%
#' @export
#' @format An R6class generator object.
tpgrp <-
  R6::R6Class(
    classname = "tpgraph",
    public = list(
      # Public Fields ----------------------------------------------------------
      nnodes = 20,
      adjacency = NULL,
      hop = "matrix",

      dim = "numeric",
      weights = "matrix",

      # Public Methods ----------------------------------------------------------
      initialize = function(topology,dim=NULL,weights=NULL){
        if(missing(topology) && is.null(self$adjacency))
          stop("to instancing the class, you must specify a topology.")

        if(!missing(topology))
          self$setTopology(topology)

        if(!is.null(dim))
          self$embedGraph(dim)

        if(!is.null(weights))
          self$weights <- weights
      },

      embedGraph = function(dim){
        # sets dimension
        if(is.null(dim)){
          self$dim <- min(sapply(1:self$nnodes,function(n)length(self$calcNeighbor(n))))-1
          message(paste("embedding dimension is automatically setted as",self$dim))
        }
        else self$dim <- dim

        self$weights <- igraph::graph.adjacency(self$adjacency) %>>%
          igraph::layout_with_kk(dim = self$dim)
      },

      setTopology = function(topology){

        # calculates adjacency matrix
        if(class(topology)=="list"){
          self$nnodes <- length(topology)
          self$adjacency <- sapply(1:self$nnodes,function(n){
            adj <- numeric(self$nnodes)
            adj[topology[[n]]] <- 1
            return(adj)
          })
        }
        else if(class(topology)=="matrix"){
          self$nnodes <- nrow(topology)
          self$adjacency <- topology
        }

        # calculates number of the hops matrix
        self$hop <- sapply(1:self$nnodes,function(center){
          h <- numeric(self$nnodes)
          i <- 1
          tip <- center
          while(!prod(h[-center])){
            tip <- which(colSums(self$adjacency[tip,,drop=F])>0)
            h[replace(logical(self$nnodes),tip,T) & h==0] <- i
            if(i==100) stop("iteration over 100 times")
            i <- i + 1
          }
          h[center] <- 0
          h
        })

      },

      calcWinner = function(x){
        if(!(class(x) %in% c("numeric","integer"))) stop("x must be numeric vector")
        if(length(x) != self$dim) stop("unmached dim and length of x")
        which.min(private$rowNorm(private$rowMinus(self$weights,x)))
      },

      calcNeighbor = function(node, neighbor.hop=1) which(self$hop[node,] <= neighbor.hop),

      calcHop = function(center,node=1:self$nnodes) self$hop[cbind(center,node)]
    ),
    private = list(
      rowNorm = function(X) sqrt(rowSums(X^2)),
      rowMinus = function(X,a) t(t(X)-a)
    )
  )
