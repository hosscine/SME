
#' Calculates coordinates on manifold.
#'
#' @param X matrix of high dimensional vectors.
#' @param A origin of manifold axis.
#' @param V direction of manifold axis.
#'
p <- function(X,A,V) rowSums((X-A) * V) / rowSums(V^2)

#' Calculates shift vector on manifold origin.
#'
#' @param n numeric of origin.
#' @param grp graph of manifold.
#'
e.inner <- function(n,grp){
  if(inherits(grp,"tpgraph"))
    i <- apply(grp$weights[grp$calcNeighbor(n,neighbor.hop=1),],2,mean) - grp$weights[n,]
  else stop(paste("unknown class",class(grp)))
  i/vnorm(i)
}

#' Calculates self-embedding model.
#'
#' @param X a matrix of data.
#' @param som a SOM learned X.
#' @param grp graph of manifold equal to SOM topology.
#'
#' @importFrom myfs matlist
#' @importFrom myfs vnorm
#'
#' @return SME object.
#' @export
#'
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

