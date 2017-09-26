#' @importFrom grDevices n2mfrow rainbow gray.colors adjustcolor
#' @importFrom graphics abline boxplot lines matlines matplot par plot polygon points
#' @importFrom utils menu combn txtProgressBar setTxtProgressBar getFromNamespace tail
#' @importFrom stats as.dist cutree hclust model.matrix predict quantile splinefun terms
#' @importFrom fda create.bspline.basis Data2fd eval.fd int2Lfd
#' @import qrcm cluster

# clusters of effects in pseudo multi-type response model: applications and theory

#' @export
clustEff <- function(Beta, p, alpha, k, ask=FALSE, k.min=1, k.max=min(10, (ncol(Beta)-1)),
                     cluster.effects=TRUE, Beta.lower=NULL, Beta.upper=NULL,
                     step=c("both", "shape", "distance"), plot=TRUE, approx.spline=FALSE, nbasis=50,
                     method="ward.D2"){

  if(!is.matrix(Beta)) Beta <- as.matrix(Beta)
  n <- nrow(Beta)
  q <- ncol(Beta)
  if(ncol(Beta) < 2) stop("At least two curves to use this procedure!")
  cut.method <- if(cluster.effects) "length" else "mindist"
  if(!is.null(Beta.lower) & !is.null(Beta.lower) & cluster.effects) cut.method <- "conf.int"
  if(cut.method == "conf.int"){
    code <- 0
    if((nrow(Beta.lower) != n) | (nrow(Beta.upper) != n)) code <- 1
    if((ncol(Beta.lower) != q) | (ncol(Beta.upper) != q)) code <- 2
    if(code != 0) stop("Dimensions of the matrices mismatched!")
  }
  if(missing(alpha)){
    alpha <- if(cluster.effects) .25 else .50
  }else{
    if(alpha < 0 | alpha > 1) stop("alpha must be in (0,1)!")
  }
  if(missing(p)) p <- 1:n/n
  if(approx.spline){
    splines <- create.bspline.basis(rangeval=range(p), nbasis=nbasis, norder=4)
    fd <- Data2fd(Beta, argvals=p, basisobj=splines)
    temp.p <- seq(fd$basis$rangeval[1], fd$basis$rangeval[2], l=max(c(501, 50*10+1)))
    Beta <- eval.fd(temp.p, fd, int2Lfd(0))
    if(cut.method == "conf.int"){
      fd.lower <- Data2fd(Beta.lower, argvals=p, basisobj=splines)
      fd.upper <- Data2fd(Beta.upper, argvals=p, basisobj=splines)
      Beta.lower <- eval.fd(temp.p, fd.lower, int2Lfd(0))
      Beta.upper <- eval.fd(temp.p, fd.upper, int2Lfd(0))
    }
    p <- temp.p
  }

  rangeBeta <- range(Beta)
  lenP <- length(p)
  METHODS <- c("ward.D", "single", "complete", "average", "mcquitty",
               "median", "centroid", "ward.D2")
  i.meth <- pmatch(method, METHODS)
  if(is.na(i.meth)) stop("invalid clustering method", paste("", method), "\n use one of this methods:", paste("", METHODS))
  if(i.meth == -1) stop("ambiguous clustering method", paste("", method))
  step <- match.arg(step)
  id <- combn(1:q, 2)
  mat <- matrix(NA, q, q)
  matDist <- matrix(NA, ncol(id), lenP)
  alphaDist <- NULL
  # BetaScale <- scale(Beta)
  distance <- vector(length=max(1, length(k.min:k.max)))
  BetaUpper <- BetaLower <- NULL

  cat("\nStep 1 (calculating alpha-percentiles): \n", sep="")
  if(step != "shape"){
    pb <- txtProgressBar(min=0, max=ncol(id), style=3)
    for(i in 1:ncol(id)){
      setTxtProgressBar(pb, i)
      s11 <- Beta[, id[1, i]]
      s22 <- Beta[, id[2, i]]
      matDist[i, ] <- abs(s11-s22)
    }
    alphaDist <- apply(matDist, 2, function(.x) quantile(.x, alpha))
  }else{
    pb <- txtProgressBar(min=0, max=1, style=3)
    setTxtProgressBar(pb, 1)
  }
  close(pb)

  cat("Step 2 (calculating dissimilarity matrix): \n", sep="")
  pb <- txtProgressBar(min=0, max=ncol(id), style=3)
  for(i in 1:ncol(id)){
    setTxtProgressBar(pb, i)
    if(step != "distance"){
      s1 <- splinefun(p, Beta[, id[1, i]])
      s2 <- splinefun(p, Beta[, id[2, i]])
      tem <- 1*I(sign(s1(p, deriv=2))*sign(s2(p, deriv=2)) == 1)
      if(step != "both") tem2 <- 0
    }

    if(step != "shape"){
      tem2 <- 1*I(matDist[i, ] <= alphaDist)
      if(step != "both") tem <- 0
    }

    mat[id[1, i], id[2, i]] <- switch(step,
                                      "both"=mean(tem & tem2),#mean(tem)*mean(tem2),#mean(tem & tem2),
                                      "distance"=mean(tem2),
                                      "shape"=mean(tem))
  }
  close(pb)

  mat <- 1-as.dist(t(mat))
  ogg <- hclust(mat, method=method)
  if(q != 2){
    if(plot) par(mfrow=c(2, 2)) else par(mfrow=c(1, 1))
    plot(ogg, main="", xlab="", ann=TRUE, sub="")
  }else{
    if(plot) par(mfrow=c(2, 1))
    k <- 2
  }

  if(missing(k)){
    if(ask){
      pick <- 1
      while(pick > 0 && pick < (q+1)){
        pick <- menu(1:q, title = "\n Select how many clusters (or 0 to exit):\n")
        if(pick > 0 && pick < (q+1)){
          break
        }
        if(pick == 0) stop("Cluster algorithm aborted!!!")
      }
    }else{
      cat("Step 3 (choosing the optimal number of clusters): \n", sep="")
      pick <- switch(cut.method,
                     "mindist"={
                       pb <- txtProgressBar(min=0, max=max(1, k.max-k.min), style=3)
                       for(pick in k.min:k.max){
                         setTxtProgressBar(pb, pick)

                         oc <- cutree(ogg, k=pick)
                         num.clust <- length(unique(oc))

                         # BetaDistMedio <- sapply(seq_len(num.clust), function(.i){
                         #   tt <- which(oc == .i)
                         #   tempcont <- NULL
                         #   if(length(tt) == 1){
                         #     tempcont <- c(tempcont, 1)
                         #   }else{
                         #     tempid <- combn(which(oc == .i), 2)
                         #     for(i in 1:ncol(tempid)) tempcont <- c(tempcont, 1-as.matrix(mat)[tempid[1, i], tempid[2, i]])
                         #     # mean(tempcont)
                         #   }
                         #   mean(tempcont)
                         # })
                         #
                         # distance[pick] <- mean(BetaDistMedio)

                         BetaMedio <- matrix(NA, nrow=lenP, ncol=num.clust)
                         BetaMedio[, 1] <- apply(as.matrix(Beta[, which(oc == 1)]), 1, mean)
                         if(num.clust > 1){
                           for(i in 2:num.clust){
                             BetaMedio[, i] <- apply(as.matrix(Beta[, which(oc == i)]), 1, mean)
                           }
                         }

                         BetaDist <- sapply(1:num.clust, function(.i) abs(Beta[, oc == .i] - BetaMedio[, .i]), simplify=FALSE)
                         idClust <- as.integer(names(which(table(oc) > 1)))
                         temp <- sort(unlist(sapply(1:length(idClust), function(.i) which(oc == idClust[.i]))))
                         BetaDistMedio <- sapply(as.integer(names(table(oc[temp]))), function(.i){
                           colM <- colMeans(BetaDist[[.i]]) # integrali
                           colM}, simplify=FALSE)
                           # 1-colM/max(colM)}, simplify=FALSE)

                         # distance[pick] <- max(unlist(lapply(BetaDistMedio, function(.x) max(.x[-which.min(.x)]))))
                         distance[pick] <- max(unlist(lapply(BetaDistMedio, function(.x) max(.x))))
                       }
                       close(pb)

                       names(distance) <- k.min:k.max
                       if(distance[1] > 1.5*distance[2]) distance <- distance[-1]

                       tempIndex2 <- which.min(diff(distance))
                       as.integer(names(tempIndex2))},
                     "length"={
                       pb <- txtProgressBar(min=0, max=1, style=3)
                       setTxtProgressBar(pb, 1)
                       distance <- diff(rev(ogg$height))
                       close(pb)

                       names(distance) <- 2:(length(distance)+1)
                       tempIndex <- which(distance < 0)
                       tempIndex2 <- which.min(distance[tempIndex])

                       as.integer(names(tempIndex2))},
                     "conf.int"={
                       pb <- txtProgressBar(min=0, max=max(1, k.max-k.min), style=3)
                       for(pick in k.min:k.max){
                         setTxtProgressBar(pb, pick)

                         oc <- cutree(ogg, k=pick)
                         num.clust <- length(unique(oc))

                         tempDist <- vector(length=num.clust)
                         for(.i in unique(oc)){
                           indCl <- (oc == .i)
                           BetaLmean <- apply(as.matrix(Beta.lower[, indCl]), 1, mean)
                           BetaRmean <- apply(as.matrix(Beta.upper[, indCl]), 1, mean)

                           xx1 <- as.matrix(BetaRmean - Beta[, indCl])
                           xx2 <- as.matrix(Beta[, indCl] - BetaLmean)

                           tempDist[.i] <- mean(colSums(apply(xx1, 2, function(.x) .x < 0))/lenP) +
                                           mean(colSums(apply(xx2, 2, function(.x) .x < 0))/lenP)
                         }

                         distance[pick] <- sum(tempDist*table(oc))/q
                        }
                        close(pb)

                        names(distance) <- k.min:k.max
                        tempIndex2 <- which.min(diff(distance))
                        as.integer(names(tempIndex2))})
    }
  }else{
    pick <- k
  }

  if(q != 2) abline(h=mean(c(max(ogg$height), rev(ogg$height), 0)[c(pick, (pick+1))]), col=2, lty=2)
  oc <- cutree(ogg, k=pick)
  num.clust <- length(unique(oc))

  Signif.interval <- NULL
  BetaMedio <- matrix(NA, nrow=lenP, ncol=num.clust, dimnames=list(p, paste0("cluster", 1:num.clust)))
  BetaMedio[, 1] <- apply(as.matrix(Beta[, which(oc == 1)]), 1, mean)
  if(!is.null(Beta.lower) & !is.null(Beta.upper)){
    Signif.interval <- BetaUpper <- BetaLower <- matrix(NA, nrow=lenP, ncol=num.clust, dimnames=list(p, paste0("cluster", 1:num.clust)))
    BetaLower[, 1] <- apply(as.matrix(Beta.lower[, which(oc == 1)]), 1, mean)
    BetaUpper[, 1] <- apply(as.matrix(Beta.upper[, which(oc == 1)]), 1, mean)
    Signif.interval[, 1] <- apply(cbind(BetaLower[, 1], BetaUpper[, 1]), 1, prod)
  }
  if(plot){
    matplot(p, Beta[, which(oc == 1)], type="l", col=1, lty=1, ylim=rangeBeta, lwd=.3, ylab="s(p)")
    lines(p, BetaMedio[, 1], col=1, lwd=1)
  }
  if(num.clust > 1){
    for(i in 2:num.clust){
      BetaMedio[, i] <- apply(as.matrix(Beta[, which(oc == i)]), 1, mean)
      if(!is.null(Beta.lower) & !is.null(Beta.upper)){
        BetaLower[, i] <- apply(as.matrix(Beta.lower[, which(oc == i)]), 1, mean)
        BetaUpper[, i] <- apply(as.matrix(Beta.upper[, which(oc == i)]), 1, mean)
        Signif.interval[, i] <- apply(cbind(BetaLower[, i], BetaUpper[, i]), 1, prod)
      }
      if(plot){
        matlines(p, Beta[, which(oc == i)], col=i, lwd=.3, lty=i)
        lines(p, BetaMedio[, i], col=i, lwd=1)
      }
    }
  }

  # BM <- apply(Beta, 1, mean)
  # BDM <- colMeans(abs(Beta - BM))

  BetaDist <- sapply(1:num.clust, function(.i) abs(Beta[, oc == .i] - BetaMedio[, .i]), simplify=FALSE)
  idClust <- as.integer(names(which(table(oc) > 1)))
  temp <- sort(unlist(sapply(1:length(idClust), function(.i) which(oc == idClust[.i]))))
  oggSil <- if(length(unique(oc[temp])) > 1) silhouette(oc[temp], dmatrix=as.matrix(mat)[temp, temp]) else NULL

  BetaDistMedio <- sapply(as.integer(names(table(oc[temp]))), function(.i){
    colM <- colMeans(BetaDist[[.i]])
    #temp.ind <- match(names(colM), names(BDM))
    colM}, simplify=FALSE)#/BDM[temp.ind]
    # 1-colM/max(colM)}, simplify=FALSE)

  DissMedio <- sapply(seq_len(num.clust), function(.i){
    tt <- which(oc == .i)
    tempcont <- NULL
    if(length(tt) == 1){
      tempcont <- c(tempcont, 1)
    }else{
      tempid <- combn(which(oc == .i), 2)
      for(i in 1:ncol(tempid)) tempcont <- c(tempcont, 1-as.matrix(mat)[tempid[1, i], tempid[2, i]])
    }
    tempcont
  }, simplify=FALSE)

  if(plot & (length(BetaDistMedio) > 0)) boxplot(DissMedio, names=as.integer(names(table(oc))))#, ylim=c(0,1))
  if(plot & (ncol(BetaMedio) > 1)) matplot(p, BetaMedio, type="l", lwd=1, ylab="s(p)")
  if(plot) par(mfrow=c(1,1))

  nms <- colnames(Beta)
  if(is.null(nms)) nms <- paste0("X", 1:q)
  colnames(Beta) <- nms
  obj <- list(call=match.call(), X=Beta, X.mean=BetaMedio, X.mean.dist=BetaDistMedio, X.lower=Beta.lower,
              X.mean.lower=BetaLower, X.upper=Beta.upper, X.mean.upper=BetaUpper, Signif.interval=Signif.interval,
              k=pick, p=p, diss.matrix=mat, X.mean.diss=DissMedio, oggSilhouette=oggSil, oggHclust=ogg, clusters=oc,
              alpha.dist=alphaDist, distance=distance, step=step, method=method, cut.method=cut.method,
              alpha=alpha)
  class(obj) <- "clustEff"

  return(obj)
}

#' @export
print.clustEff <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  tabClust <- as.integer(table(x$clusters))
  cat("Clustering Effects algorithm with ", x$k, " clusters of sizes ",
      paste(tabClust, collapse = ", "), "\n", sep = "")
  cat("\n")

  invisible(x)
}

#' @export
summary.clustEff <- function(object, ...){
  nClust <- as.integer(table(object$clusters))

  avClust <- vector(length=object$k)
  avClust[which(nClust > 1)] <- unlist(lapply(object$X.mean.dist, mean))
  names(avClust) <- paste0("Cluster", seq(object$k))

  avDiss <- vector(length=object$k)
  avDiss <- unlist(lapply(object$X.mean.diss, mean))
  names(avDiss) <- paste0("Cluster", seq(object$k))

  if(!is.null(object$oggSilhouette)){
    avSilhouette <- rep(1, object$k)
    avSilhouette[which(nClust > 1)] <- summary(object$oggSilhouette)$clus.avg.widths
    names(avSilhouette) <- paste0("Cluster", seq(object$k))
  }else{
    avSilhouette <- NULL
  }

  obj <- list(call=object$call, k=object$k, n=nrow(object$X), p=ncol(object$X), step=object$step, alpha=object$alpha,
              method=object$method, cut.method=object$cut.method, tabClust=table(object$clusters),
              avClust=avClust, avSilhouette=avSilhouette, avDiss=avDiss)

  class(obj) <- "summary.clustEff"
  return(obj)
}

#' @export
print.summary.clustEff <- function(x, digits=max(3L, getOption("digits") - 3L), ...){

  cat("\n######################", "\n\n")
  cat("Selected n. of clusters:", format(x$k, digits=digits), "\n")
  cat("n. of observations:", x$n, "\n")
  cat("n. of curves:", x$p, "\n")
  cat("Step selected:", x$step, "\n")
  cat("alpha-percentile selected:", x$alpha, "\n")
  cat("Cut method selected:", x$cut.method, "\n")
  cat("Clustering method selected:", x$method, "\n\n")
  cat("######################", "\n\n")

  cat("Clustering table:")
  print(x$tabClust, ...)
  cat("\n######################", "\n\n")

  cat("Average cluster distance:\n")
  print(round(x$avClust, digits), ...)
  cat("\n######################", "\n\n")

  cat("Average cluster dissimilarity:\n")
  print(round(x$avDiss, digits), ...)
  cat("\n######################", "\n\n")

  if(!is.null(x$avSilhouette)){
    cat("Individual silhouette widths:\n")
    print(round(x$avSilhouette, digits), ...)
    cat("\n######################", "\n\n")
  }

  invisible(x)
}

#' @export
plot.clustEff <- function(x, xvar=c("clusters", "dendrogram", "boxplot"), which, add=FALSE,
                          all=FALSE, polygon=TRUE, dissimilarity=TRUE, ...){
  xvar = match.arg(xvar)
  if(!missing(which)){
    if(any(which <= 0) | any(which > x$k)){
      stop("Which values in 1-", x$k)
    }
  }

  L <- list(...)

  switch(xvar,
         "dendrogram"={
           if(is.null(L$main)) L$main <- "Dendrogram"
           if(is.null(L$col)) L$col <- "red"
           if(is.null(L$xlab)) L$xlab <- ""
           if(is.null(L$lty)) L$lty <- 2

           # par(mfrow=c(1,1))
           plot(x$oggHclust, main=L$main, xlab=L$xlab, ann=TRUE, sub="")
           abline(h=mean(c(max(x$oggHclust$height), rev(x$oggHclust$height), 0)[c(x$k, (x$k+1))]), col=L$col, lty=L$lty)
         },
         "clusters"={
           if(missing(which)) which <- seq(x$k)
           # if(!missing(which) & !add & (length(which) != 1)) par(mfrow=n2mfrow(length(which)))# else par(mfrow=c(1,1))
           if(all & !add) par(mfrow=n2mfrow(length(which)))
           yRange <- range(x$X[,  x$clusters %in% which])
           if(!is.null(x$X.lower) & !is.null(x$X.upper)){
              temp <- range(cbind(x$X.mean.lower, x$X.mean.upper))
              yRange <- range(c(temp, yRange))
           }

           tabClust <- table(x$clusters)

           if(is.null(L$xlab)) L$xlab <- "p"
           if(is.null(L$ylab)) L$ylab <- "s(p)"
           if(is.null(L$main)) L$main <- ""
           if(is.null(L$lwd)) L$lwd <- c(1, 1.5)
           if(length(L$lwd) == 1) L$lwd <- rep(L$lwd, 2)
           if(is.null(L$type)) L$type <- "l"
           if(is.null(L$ylim)) L$ylim <- yRange
           if(is.null(L$lty)){L$lty <- 1}
           # if(is.null(L$legend)) L$legend <- "topright"
           if(is.null(L$col)) L$col <- gray.colors(length(which))
           if(length(L$col) == 1) L$col <- rep(L$col, length(which))

           if(length(which) == 1){
             tempInd <- x$clusters == which
             matplot(x$p, x$X[, tempInd], xlab=L$xlab, ylab=L$ylab, type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=1, lty=L$lty)
             lines(x$p, x$X.mean[, which], lwd=L$lwd[2], col=2)
             # if(!is.null(x$X.lower) & !is.null(x$X.upper)){
             #   lines(x$p, x$X.mean.lower[, which], lwd=L$lwd[2], col=2, lty=2)
             #   lines(x$p, x$X.mean.upper[, which], lwd=L$lwd[2], col=2, lty=2)
             # }
             if(!is.null(x$X.lower) & !is.null(x$X.upper)){
               if(polygon){
                 yy <- c(x$X.mean.lower[, which], tail(x$X.mean.upper[, which], 1), rev(x$X.mean.upper[, which]),
                         x$X.mean.lower[1, which])
                 xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
                 polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
               }else{
                 points(x$p, x$X.mean.lower[, which], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
                 points(x$p, x$X.mean.upper[, which], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
               }
             }
             abline(h=0, lty=3, col=1)
             # legend(L$legend, legend=colnames(x$X)[tempInd], bty="n", lty=seq(tempInd), col=1)
           }else{
             matplot(x$p, x$X[, x$clusters == which[1]], xlab=L$xlab, ylab=L$ylab, type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty)
             lines(x$p, x$X.mean[, which[1]], lwd=L$lwd[2], col=2)
             if(!is.null(x$X.lower) & !is.null(x$X.upper)){
              if(polygon){
                yy <- c(x$X.mean.lower[, which[1]], tail(x$X.mean.upper[, which[1]], 1), rev(x$X.mean.upper[, which[1]]), x$X.mean.lower[1, which[1]])
                xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
                polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
              }else{
                points(x$p, x$X.mean.lower[, which[1]], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
                points(x$p, x$X.mean.upper[, which[1]], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
              }
             }
             abline(h=0, lty=3, col=1)

             # if(!is.null(x$X.lower) & !is.null(x$X.upper)){
             #   lines(x$p, x$X.mean.lower[, which[1]], lwd=L$lwd[2], col=2, lty=2)
             #   lines(x$p, x$X.mean.upper[, which[1]], lwd=L$lwd[2], col=2, lty=2)
             # }


             for(i in 2:length(which)){
               if(!add){
                 matplot(x$p, x$X[, x$clusters == which[i]], xlab=L$xlab, ylab=L$ylab, type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[i], lty=L$lty)
               }else{
                 matlines(x$p, x$X[, x$clusters == which[i]], lwd=L$lwd[1], col=L$col[i], lty=L$lty)
               }
               lines(x$p, x$X.mean[, which[i]], lwd=L$lwd[2], col=2)
               if(!is.null(x$X.lower) & !is.null(x$X.upper)){
                 if(polygon){
                   yy <- c(x$X.mean.lower[, which[i]], tail(x$X.mean.upper[, which[i]], 1), rev(x$X.mean.upper[, which[i]]), x$X.mean.lower[1, which[i]])
                   xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
                   polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
                 }else{
                   points(x$p, x$X.mean.lower[, which[i]], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
                   points(x$p, x$X.mean.upper[, which[i]], lty = 2, lwd = L$lwd[2], type = "l", col = 2)
                 }
               }
               # if(!is.null(x$X.lower) & !is.null(x$X.upper)){
               #   lines(x$p, x$X.mean.lower[, which[i]], lwd=L$lwd[2], col=2, lty=2)
               #   lines(x$p, x$X.mean.upper[, which[i]], lwd=L$lwd[2], col=2, lty=2)
               # }
               abline(h=0, lty=3, col=1)
             }
           }

           #par(mfrow=c(1,1))
           },
         "boxplot"={
           tabClust <- table(x$clusters)
           X <- if(dissimilarity) x$X.mean.diss else x$X.mean.dist

           if(is.null(L$main)) L$main <- "Average cluster distance"
           if(is.null(L$labels)) L$labels <- as.integer(names(tabClust))[tabClust > 1]
           if(is.null(L$ylab)) L$ylab <- "Distance"
           if(is.null(L$ylim)) L$ylim <- c(0, 1)
           if(!dissimilarity) L$ylim <- c(0, max(unlist(X)))

           # par(mfrow=c(1,1))
           if(length(X) > 0){
             boxplot(X, names=L$labels, ylim=L$ylim, main=L$main, ylab=L$ylab)
           }
         })
}

#' @export
extract.object <- function(Y, X, intercept=TRUE, formula.p=~slp(p, 3), s, object, p, which){
  if(missing(p)) p <- seq(.01, .99, .01)

  if(!missing(object)){
    if(!inherits(object, "iqr")){
      stop("Wrong class object!")
    }else{
      if(is.data.frame(object$mf)){
        X <- as.matrix(object$mf[, -1])
      }else{
        X <- object$mf[[2]]
      }

      n <- nrow(X)
      q <- ncol(X)
      intercept <- attr(attr(object$mf, "terms"), "intercept")
      if(missing(which)) which <- 1:q
      if(intercept) which <- which + 1
      labels <- colnames(X)

      tempX <- tempXl <- tempXr <- list()
      index <- 0
      for(i in which){
        index <- index + 1
        predObj <- predict(object, type="beta", p=p)
        tempX[[index]] <- predObj[[i]][, 2]
        tempXl[[index]] <- predObj[[i]][, 4]
        tempXr[[index]] <- predObj[[i]][, 5]
      }

      names(tempX) <- names(tempXl) <- names(tempXr) <- labels
    }
  }else{
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    qY <- ncol(Y)
    qX <- ncol(X) + intercept
    labels <- colnames(X)
    if(is.null(labels)){
      labels <- paste0(X, 1:(qX-1))
      colnames(X) <- labels
    }
    labels <- if(intercept) c("(Intercept)", labels) else labels
    n <- nrow(X)
    if(n != nrow(Y)) stop("Dimension  mismatched!")
    is.slp <- getFromNamespace("is.slp", "qrcm")
    if(use.slp <- is.slp(formula.p)){
      k <- attr(use.slp, "k")
      intB <- attr(use.slp, "intB")
      termlabelsB <- paste("slp", 1:k, sep = "")
      k <- k + intB
      coefnamesB <- (if (intB) c("(Intercept)", termlabelsB) else termlabelsB)
    }else{
      B <- model.matrix(formula.p, data = data.frame(p=p))
      k <- ncol(B)
      termlabelsB <- attr(terms(formula.p), "term.labels")
      coefnamesB <- colnames(B)
    }
    if(missing(s)) s <- matrix(1, nrow=qX, ncol=k)
    colnames(s) <- coefnamesB
    rownames(s) <- labels
    if(missing(which)) which <- if(intercept) 1:(qX-1) else 1:qX
    if(intercept) which <- which + 1

    tempX <- tempXl <- tempXr <- list()
    index <- 0
    for(i in which){
      index <- index + 1
      tempXr[[index]] <- tempXl[[index]] <- tempX[[index]] <- matrix(NA, length(p), qY)
    }

    for(i in 1:qY){
      object <- if(intercept) iqr(Y[, i] ~ X, formula.p=formula.p, s=s) else iqr(Y[, i] ~ -1 + X, formula.p=formula.p, s=s)
      index <- 0
      for(j in which){
        index <- index + 1
        predObj <- predict(object, type="beta", p=p)
        tempX[[index]][, i] <- predObj[[j]][, 2]
        tempXl[[index]][, i] <- predObj[[j]][, 4]
        tempXr[[index]][, i] <- predObj[[j]][, 5]
      }
    }

    names(tempX) <- names(tempXl) <- names(tempXr) <- if(intercept) labels[-1] else labels
  }

  return(list(p=p, X=tempX, Xl=tempXl, Xr=tempXr))
}
