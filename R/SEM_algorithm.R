
#' calculate coordinates on manifold
#' @param X Matrix of high dimensional vectors.
#' @param A Origin of manifold axis.
#' @param V Direction of manifold axis.
p <- function(X,A,V) rowSums((X-A) * V) / rowSums(V^2)

#' calculate shift vector on manifold origin
#' @param n Numeric of origin.
#' @param grp Graph of manifold.
e.inner <- function(n,grp){
  if(inherits(grp,"tpgraph"))
    i <- apply(grp$weights[grp$calcNeighbor(n,neighbor.hop=1),],2,mean) - grp$weights[n,]
  else stop(paste("unknown class",class(grp)))
  i/vnorm(i)
}

#' calculate self-embedding model
#' @export
#' @param X A matrix of data.
#' @param som A SOM learned X.
#' @param grp Graph of manifold equal to som topology.
#' @return SMM object
#' @import pipeR
#' @examples library(TDA)
#'
#' sphere4d <- sphereUnif(n = 500,d = 3)
#' spsom <- setSom(4,topo.mode = "sphere")
#' spgrp <- setGraph(spsom)
#'
#' spsem <- sem(sphere4d,spsom,spgrp)
sem <- function(X,som,grp){
  if(!is.matrix(X)){X <- rbind(X,X);single <- T}
  else single <- F
  cl <- apply(X,1,som$calcWinner)
  emdim <- grp$dim - 1

  # select node using for stretch boudary simplex
  cl.nei <- lapply(1:grp$nnodes,grp$calcNeighbor)
  nei.max <- max(sapply(cl.nei,length))
  cl.nei <- matlist(lapply(cl.nei,function(nei){length(nei)<-nei.max;nei[2:nei.max]}))
  nareplace <- is.na(cl.nei) * 1:som$nnodes
  cl.nei[is.na(cl.nei)] <- 0
  cl.nei <- cl.nei + nareplace
  nei.max <- nei.max - 1

  # edge position
  A <- som$weights[cl,]
  a <- grp$weights[cl,]
  pV <- matrix(0,nrow(A),ncol(A))
  pv <- matrix(0,nrow(a),ncol(a))

  allpar <- sapply(1:nei.max,function(n){
    V <- som$weights[cl.nei[cbind(cl,n)],] - A
    p(X,A,V)
  })
  allpar[is.nan(allpar)] <- 0
  abspar <- abs(allpar)

  par <- matrix(0,nrow(X),emdim)
  simplex.node <- matrix(nrow = nrow(X),ncol = emdim)
  for(i in 1:emdim){
    parind <- cbind(1:nrow(allpar),max.col(abspar))
    simplex.node[,i] <- cl.nei[cbind(cl,parind[,2])]
    par[,i] <- allpar[parind]
    abspar[parind] <- -Inf
  }
  for(n in 1:emdim){
    pV <- pV + par[,n]*(som$weights[simplex.node[,n],] - A)
    pv <- pv + par[,n]*(grp$weights[simplex.node[,n],] - a)
  }

  Inners <- t(sapply(1:som$nnodes,e.inner,som))
  inners <- t(sapply(1:grp$nnodes,e.inner,grp))

  # direction of egde to x and distance of edge to x on the innner vector
  direction <- X - A - pV
  d <- rowSums(Inners[cl,]*direction)
  Y <- a + pv + inners[cl,]*d

  if(single){
    X <- X[1,]; Y <- Y[1,]; cl <- cl[1]; simplex.node <- simplex.node[1,];
    p <- par[1,]; d <- d[1]}

  semR6$new(X = X, Y = Y, cl = cl, delta = simplex.node, p = par,d = d,
            som = som, grp = grp)
}


#' Deforms listed list (a list nested by list) as matrix
#'
#' @param lst listed list.
#' @param nrow nrow of returned matrix.
#' @param ncol ncol of returned matrix.
#' @param rname rownames of returned matrix.
#' @param cname colnames of returned matrix.
#'
#' @return matrix deformed from listed list.
matlist = function(lst,nrow=length(lst),ncol=length(lst[[nrow]]),
                   rname=names(lst),cname=names(lst[[nrow]])){
  if (length(lst)==0) stop("input list is empty")
  lst <- lapply(lst,function(l){length(l)<-ncol;l})
  mat <- matrix(unlist(lst), nrow, ncol,byrow = T)
  rownames(mat) <- rname
  colnames(mat) <- cname
  return(mat)
}

#' calculates vector's norm
#'
#' @param x input vector.
#'
#' @return norm of \code{x}.
vnorm <- function(x) sqrt(sum(x^2))
