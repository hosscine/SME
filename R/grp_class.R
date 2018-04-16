#' R6 generator object of topological graph.
#'
#' \code{tpgrp$new(topology, dim = NULL, weights = NULL)}
#'
#' @importFrom R6 R6Class
#' @importFrom magrittr %>%
#' @importFrom rgl plot3d
#' @importFrom myfs overwriteEllipsis
#' @importFrom myfs rowNorm
#' @importFrom myfs rowMinus
#'
#' @export
#' @docType class
#' @format An R6class generator object.
#'
tpgrp <-
  R6::R6Class(
    classname = "tpgraph",
    public = list(
      # Public Fields ----------------------------------------------------------
      nnodes = 20,
      adjacency = NULL,
      hop = "matrix",

      dim = "numeric",
      weights = NULL,

      # Public Methods ----------------------------------------------------------
      initialize = function(topology, dim = NULL, weights = NULL){
        if(missing(topology) && is.null(self$adjacency))
          stop("to instancing the class, you must specify a topology.")

        if(!missing(topology))
          self$setTopology(topology)

        if(!is.null(dim) || is.null(weights))
          self$weights <- self$embedGraph(dim, setdim = T)

        if(!is.null(weights))
          self$weights <- weights
      },

      embedGraph = function(dim = self$dim, setdim = F){
        if(is.null(dim)){
          emdim <- min(sapply(1:self$nnodes, function(n)length(self$calcNeighbor(n))), 4) - 1
          message(paste("embedding dimension is automatically set as",self$dim))
        }
        else emdim <- dim

        if(setdim) self$dim <- emdim

        return(igraph::graph.adjacency(self$adjacency) %>%
                 igraph::layout_with_kk(dim = emdim))
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
        else if(class(topology) == "matrix"){
          self$nnodes <- nrow(topology)
          self$adjacency <- topology
        }

        # calculates number of the hops matrix
        self$hop <- sapply(1:self$nnodes,function(center){
          h <- numeric(self$nnodes)
          i <- 1
          tip <- center
          while(!prod(h[-center])){
            tip <- which(colSums(self$adjacency[tip,, drop = F]) > 0)
            h[replace(logical(self$nnodes), tip, T) & h == 0] <- i
            if(i == 100) stop("iteration over 100 times")
            i <- i + 1
          }
          h[center] <- 0
          h
        })

      },

      calcWinner = function(x){
        if(!(class(x) %in% c("numeric", "integer"))) stop("x must be numeric vector")
        if(length(x) != self$dim) stop("dimension of weights and x are differ")
        which.min(rowNorm(rowMinus(self$weights, x)))
      },

      calcNeighbor = function(node, neighbor.hop=1) which(self$hop[node,] <= neighbor.hop),

      calcHop = function(center,node=1:self$nnodes) self$hop[cbind(center,node)],

      plot = function(..., dim = self$dim, reemb = F){
        if(reemb) x <- self$embedGraph(dim)
        else x <- self$weights

        if(dim == 2){
          elp <- overwriteEllipsis(..., xlab = "", ylab = "", x = x)
          do.call(plot,elp)
          private$drawNodeEdges(dim = dim, reemb = reemb)
        }
        else if(dim == 3){
          elp <- overwriteEllipsis(..., xlab = "", ylab = "", zlab = "", x = x)
          do.call(plot3d,elp)
          private$drawNodeEdges(dim = dim, reemb = reemb)
        }
        else{
          grp <- igraph::graph.adjacency(self$adjacency)
          plot(grp, mode = "undirected", ...)
        }
      }
    ),

    private = list(
      drawNodeEdges = function(dim = self$dim, reemb = F){
        if(reemb) X <- self$embedGraph(dim)
        else if(dim != self$dim){
          X <- self$embedGraph(dim)
          cat("graph is auto re-embedded due to given dimenson",dim,"and self dimension",self$dim,"\n")
        }
        else X <- self$weights

        # drawing edges
        if(dim==2){
          for(i in 1:self$nnodes){
            nei <- self$calcNeighbor(i, neighbor.hop = 1)
            graphics::segments(X[rep(i, length(nei)), 1], X[rep(i, length(nei)), 2], X[nei, 1], X[nei, 2])
          }
          graphics::text(X, labels=1:self$nnodes, pos = 1)
        }
        else if(dim==3){
          for(i in 1:self$nnodes){
            nei <- self$calcNeighbor(i, neighbor.hop = 1)
            rgl::segments3d(X[replace(rep(i ,2 * length(nei)), 2 * 1:length(nei), nei),])
          }
          rgl::text3d(X, texts = 1:self$nnodes, adj = 1)
        }
        else stop("edges can not be drawned for dimension ", dim)
      }
    )
  )
