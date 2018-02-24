# s <- sem(data,som,graph)
# Y = s$yw[s$cl,] + s$pv + s$ys[s$cl,]*s$d
# y = sy + epsilon * theta()

#' R6 class self-embedding model
#'
#' @docType class
#' @export
semR6 <- R6::R6Class(
  classname = "sem",

  public = list(
    X = "matrix",
    Y = "matrix",
    xdim = "numeric",
    ydim = "numeric",

    cl = "numeric",

    delta = "matrix",
    ndelta = "numeric",
    p = "matrix",
    d = "numeric",

    nnodes = "numeric",
    som = "tpsom",
    grp = "tpgraph",

    initialize = function(X,Y,cl,delta,p,d,som,grp){
      self$X <- X
      self$Y <- Y
      self$xdim <- ncol(X)
      self$ydim <- ncol(Y)

      self$cl <- cl

      self$delta <- delta
      self$ndelta <- ncol(delta)
      self$p <- p
      self$d <- d

      self$nnodes <- som$nnodes
      self$som <- som
      self$grp <- grp
    },

    pvalue = function(x){
      if(is.matrix(x)){
        cl <- apply(x,1,self$som$calcWinner)
        print(cl)
        rowSums((x-self$xw[cl,]) * self$xv[cl,]) / rowSums(self$xv[cl,]^2)
      }
      else{
        cl <- apply(x,1,self$som$calcWinner)
        sum((x-self$xw[cl,]) * self$xv[cl,]) / sum(self$xv[cl,]^2)
      }
    }

  ),

  active = list(

    xv = function() private$v.direction(self$xdim,self$xw),
    yv = function() private$v.direction(self$ydim,self$yw),

    xpv = function() private$pTimesv(self$xv,self$xdim),
    ypv = function() private$pTimesv(self$yv,self$ydim),

    xs = function() t(sapply(1:self$nnodes,private$e.inner,self$som)),
    ys = function() t(sapply(1:self$nnodes,private$e.inner,self$grp)),

    xw = function() self$som$weights,
    yw = function() self$grp$weights

  ),

  private = list(

    v.direction = function(d,w){
      ret <- array(0,dim = c(nrow(self$X),d,self$ndelta))
      for(k in 1:self$ndelta)
        ret[,,k] <- w[self$delta[,k],] - w[self$cl,]
      return(ret)
    },

    pTimesv = function(v,d){
      # suppressMessages(require(pipeR))
      sapply(1:self$ndelta,FUN = function(n)self$p[,n] * v[,,n]) %>>%
        array(dim = c(nrow(self$X),d,self$ndelta)) %>>%
        apply(MARGIN = c(1,2),FUN = sum)
    },

    e.inner = function(n,grp){
      i <- apply(grp$weights[grp$calcNeighbor(n,neighbor.hop=1),],2,mean) - grp$weights[n,]
      i/vnorm(i)
    }

  )
)
