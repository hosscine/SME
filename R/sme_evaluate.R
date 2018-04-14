#' Compare the embedding of sme with the conventional methods.
#'
#' @param sme self-manifolding embedding object.
#' @param k number of neighbors for Trustworthiness and Continuity.
#' @param embk number of neighbors for embedding of LLE and Isomap.
#' @param plot.only if \code{T}, doesn't calculate evaluation.
#' @param plot.fig if \code{T}, plot embedding result.
#' @param col color of pca.
#' @param return.embs if\code{T}, returns embedding results with listed.
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom graphics par plot
#' @importFrom stats cmdscale dist
#' @importFrom utils capture.output
#'
#' @return a matrix contains evaluation.
#' @export
#'
evaluateSme <- function(sme, k = 12, embk = 10, return.embs = F,
                        plot.only = F, plot.fig = T, col = NULL){

  # data ready
  X <- sme$X
  som <- sme$som
  dist.origin <- dist(X)
  cl <- apply(X,1,som$calcWinner)
  if(missing(col)) col <- grDevices::rainbow(som$nnodes)[apply(X, 1, som$calcWinner)]
  emdim <- min(sapply(1:som$nnodes, function(n)
    length(som$calcNeighbor(n, neighbor.hop = 1)))) - 1

  # calculates embeddings
  method.names <- c("SME","PCA","MDS","Isomap","LLE")
  capture.output(emb <- list("sem" = sme$Y,
                             "pca" = prcomp(X)[[5]][, 1:emdim],
                             "mds" = cmdscale(dist.origin, k = emdim),
                             "iso" = vegan::isomap(dist.origin, ndim = emdim, k=embk)[[1]],
                             "lle" = lle::lle(X,m = emdim, k = embk)$Y))

  # plot
  if(plot.fig){
    par(mfrow = c(2, 3))
    if(emdim == 2) for(i in 1:length(emb))
      plot(emb[[i]],col=col, xlab = "", ylab = "", main = method.names[i], cex=2)
    else if(emdim == 3) for(i in 1:length(emb))
      scatterplot3d(emb[[i]], color = col, main = method.names[i])
  }
  if(plot.only) return("only plot.")

  # calculates dist of embedding results
  dst <- lapply(emb,dist)

  # embedding measure 1: distance error
  # lower is better
  em1 <- sapply(dst,function(x)sum(abs(dist.origin-x)))

  # deforms dist to distance matrix
  dist.origin <- as.matrix(dist.origin)
  dst <- lapply(dst,as.matrix)

  # calculates rank of closeness of each data
  rnk.origin <- t(apply(dist.origin,1,order))
  rnk <- lapply(dst,function(X)t(apply(X,1,order)))

  # initialize calculation of evaluation
  #[Neighborhood Preservation in Nonlinear Projection Methods: An Experimental Study]
  NN1error <- numeric(length(rnk))
  r.minus.k <- numeric(length(rnk))
  rhat.minus.k <- numeric(length(rnk))

  # loop and loops
  # O(n*methods)
  for(i in 1:nrow(X)){
    NN.origin <- rnk.origin[i, 1:k]

    for(j in 1:length(rnk)){
      rk <- rnk[[j]]
      NN <- rk[i, 1:k]
      Uk <- setdiff(NN,NN.origin)
      Vk <- setdiff(NN.origin,NN)

      if(which(rnk.origin[i,] == 1) != which(rk[i,]==1)) NN1error[j] <- NN1error[j] + 1
      if(length(Uk) != 0) r.minus.k[j] <- r.minus.k[j] + sum(rnk.origin[i,Uk] - k)
      if(length(Vk) != 0) rhat.minus.k[j] <- rhat.minus.k[j] + sum(rk[i,Vk] - k)
    }
  }

  # embedding measure 2: Generalization errors of 1-NN classifiers
  # smaller is better
  NN1 <- NN1error/nrow(X)

  # embedding measure3: Trustworthiness T(k) と Continuity C(k)
  # larger is better
  Tk <- 1 - 2 / (nrow(X) * k * (2 * nrow(X) - 3 * k - 1)) * r.minus.k
  Ck <- 1 - 2 / (nrow(X) * k * (2 * nrow(X) - 3 * k - 1)) * rhat.minus.k

  eval <- rbind(em1, NN1, Tk, Ck)
  colnames(eval) <- method.names
  if(return.embs) attr(eval, "embs") <- emb
  eval
}
