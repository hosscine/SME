
#' Optimize graph to preserve the distance between previous and after the embedding.
#'
#' @param X target data.
#' @param som learned som.
#' @param grp target graph.
#'
#' @export
#'
grpOptimize <- function(X, som, grp){
  o <- stats::optimize(grpWeighting, interval = c(0, 10),
                X, dist(X), som, grp)
  refdist <- grpWeighting(1, X, dist(X), som, grp)
  grp$weights <- grp$weights * o$minimum
  cat("metric of graph embedded space is optimized.\n")
  cat("distance error reduced", refdist, "to", o$objective, "  improved", refdist - o$objective)
}


# GRPのweightsをt倍してSMEを計算したときの距離の誤差を返す
#' Objective function for grpOptimize
#'
#' @param t parameter.
#' @param X target data.
#' @param X.dist dist of target data.
#' @param som learned som.
#' @param grp target graph.
#'
#' @importFrom assertthat assert_that
#'
#' @return error of distance due to the embedding with multiplying weights of the graph by t.
#'
grpWeighting <- function(t, X, X.dist, som, grp){
  assert_that(is.matrix(X))
  assert_that(class(X.dist) == "dist")
  assert_that(assertthat::is.number(t))

  grp <- grp$clone()
  grp$weights <- grp$weights * t
  sme <- sem(X, som, grp)
  Y.dist <- dist(sme$Y)

  sum(abs(Y.dist - X.dist))
}
