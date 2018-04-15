#' Calculates evaluate measures of the embedding.
#'
#' @param Y after embedding data or dist.
#' @param X pre embedding data or dist.
#' @param k number of neighbors for evaluate measures.
#'
#' @return evaluate measures.
#' @export
#'
evaluateEmbedding <- function(Y, X, k = 12){

  # dist ready
  if(class(X) == "dist") dst.X <- as.matrix(X)
  else if(is.matrix(X) && nrow(X) == ncol(X)) dst.X <- X
  else if(is.matrix(X)) dst.X <- as.matrix(dist(X))
  else stop("X must be dist, dist matrix or data matrix")

  if(class(Y) == "dist") dst.Y <- as.matrix(Y)
  else if(is.matrix(Y) && nrow(Y) == ncol(Y)) dst.Y <- Y
  else if(is.matrix(Y)) dst.Y <- as.matrix(dist(Y))
  else stop("Y must be dist, dist matrix or data matrix")

  # embedding measure 1: distance error
  # lower is better
  dist.diff <- sum(abs(dst.Y - dst.X))

  # calculates rank of closeness of each data
  rnk.X <- t(apply(dst.X, 1, order))
  rnk.Y <- t(apply(dst.Y, 1, order))

  # initialize calculation of evaluation
  #[Neighborhood Preservation in Nonlinear Projection Methods: An Experimental Study]
  NN1error <- 0
  r.minus.k <- 0
  rhat.minus.k <- 0

  # loop and loops
  # O(n*methods)
  for(i in 1:nrow(dst.X)){
    NN.X <- rnk.X[i, 1:k]
    NN.Y <- rnk.Y[i, 1:k]
    Uk <- setdiff(NN.Y, NN.X)
    Vk <- setdiff(NN.X, NN.Y)

    if(which(rnk.X[i,] == 1) != which(rnk.Y[i,] == 1)) NN1error <- NN1error + 1
    if(length(Uk) != 0) r.minus.k <- r.minus.k + sum(rnk.X[i, Uk] - k)
    if(length(Vk) != 0) rhat.minus.k <- rhat.minus.k + sum(rnk.Y[i, Vk] - k)
  }

  NN1 <- NN1error / nrow(dst.X)

  # embedding measure 2: Generalization errors of 1-NN classifiers
  # smaller is better
  # embedding measure 3: Trustworthiness T(k) と Continuity C(k)
  # larger is better
  Tk <- 1 - 2 / (nrow(dst.X) * k * (2 * nrow(dst.X) - 3 * k - 1)) * r.minus.k
  Ck <- 1 - 2 / (nrow(dst.X) * k * (2 * nrow(dst.X) - 3 * k - 1)) * rhat.minus.k

  return(cbind(dist.diff, NN1, Tk, Ck))
}
