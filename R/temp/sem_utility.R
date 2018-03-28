pca.graph <- function(X,som,...){
  require(rgl)
  pca <- prcomp(rbind(X,som$weights))[[5]]
  
  args <- list(...)
  if(is.null(args$size)) args$size <- 5
  if(is.null(args$col)) col <- rainbow(som$nnodes)[apply(X,1,som$calcWinner)]
  if(is.null(args$xlab)) args$xlab <- ""
  if(is.null(args$ylab)) args$ylab <- ""
  if(is.null(args$zlab)) args$zlab <- ""
  args$x <- pca[1:nrow(X),]
  # args$x <- pca
  
  do.call(plot3d, args)
  # plot3d(pca[1:nrow(X),],col=rainbow(som$nnodes)[apply(X,1,som$calcWinner)],
         # xlab="",ylab="",zlab="",...)
  points3d(pca[(nrow(X)+1):nrow(pca),],size=5)
  for(i in 1:som$nnodes){
    nei <- som$calcNeighbor(i,neighbor.hop = 1)
    segments3d(pca[replace(rep(i,2*length(nei)),2*1:length(nei),nei)+nrow(X),])
  }
  text3d(pca[(nrow(X)+1):nrow(pca),],texts = 1:som$nnodes,adj = 1)
}

grpOptimize <- function(X, som, grp){
  o <- optimize(grpWeighting, interval = c(0,10),
                X, dist(X), som, grp)
  refdist <- grpWeighting(1, X, dist(X), som, grp)
  grp$weights <- grp$weights * o$minimum
  cat("metric of graph embedded space is optimized.\n")
  cat("distance error reduced", refdist, "to", o$objective, "  improved", refdist - o$objective)
}

# GRPのweightsをt倍してSMEを計算したときの距離の誤差を返す
#' Title
#'
#' @param t 
#' @param X 
#' @param X.dist 
#' @param som 
#' @param grp 
#'
#' @return
#' @export
#'
#' @examples
#' optimize(grpWeighting,interval = c(0,10),
#'  　　　　ahiru.nz.mat,dist(ahiru.nz.mat),ahiru.nz.som,ahiru.nz.grp)
grpWeighting <- function(t, X, X.dist, som, grp){
  assert_that(is.matrix(X))
  assert_that(class(X.dist)=="dist")
  assert_that(is.number(t))
  
  grp <- grp$clone()
  grp$weights <- grp$weights * t
  sme <- sem(X,som,grp)
  Y.dist <- dist(sme$Y)
  
  sum(abs(Y.dist - X.dist))
}
