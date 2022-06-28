#' @importFrom grDevices n2mfrow rainbow gray.colors adjustcolor gray
#' @importFrom graphics abline boxplot lines matlines matplot par plot polygon points axis mtext
#' @importFrom utils menu combn txtProgressBar setTxtProgressBar getFromNamespace tail
#' @importFrom stats as.dist cutree hclust model.matrix predict quantile splinefun terms dist cor
#' @importFrom fda create.bspline.basis Data2fd eval.fd int2Lfd pca.fd
#' @import qrcm cluster ggpubr ggplot2

# clusters of effects in pseudo multi-type response model: applications and theory

#' @export
clustEff <- function(Beta, Beta.lower = NULL, Beta.upper = NULL,
                     k = c(2, min(5, (ncol(Beta)-1))), ask = FALSE,
                     diss.mat, alpha = .5,
                     step = c("both", "shape", "distance"),
                     cut.method = c("mindist", "length", "conf.int"),
                     method = "ward.D2", approx.spline = FALSE, nbasis = 50,
                     conf.level = 0.9, stand = FALSE, plot = TRUE, trace = TRUE){

  if(!is.matrix(Beta)) Beta <- as.matrix(Beta)
  if(stand) Beta <- scale(Beta)
  n <- nrow(Beta); q <- ncol(Beta)
  if(q < 2) stop("At least two curves to use this procedure!")
  p <- 1:n / n; lenP <- length(p)
  if(approx.spline){
    splines <- create.bspline.basis(rangeval=range(p), nbasis=nbasis, norder=3)
    fd <- Data2fd(Beta, argvals=p, basisobj=splines)
    temp.p <- seq(fd$basis$rangeval[1], fd$basis$rangeval[2], l=max(c(501, 50*10+1)))
    Beta <- eval.fd(temp.p, fd, int2Lfd(0))
    if(!is.null(Beta.lower) & !is.null(Beta.upper)){
      fd.lower <- Data2fd(Beta.lower, argvals=p, basisobj=splines)
      fd.upper <- Data2fd(Beta.upper, argvals=p, basisobj=splines)
      Beta.lower <- eval.fd(temp.p, fd.lower, int2Lfd(0))
      Beta.upper <- eval.fd(temp.p, fd.upper, int2Lfd(0))
    }
    p <- temp.p
  }

  nms <- colnames(Beta)
  if(is.null(nms)) nms <- paste0("X", 1:q)
  colnames(Beta) <- nms

  if(length(k) > 2) stop("The length of k have to be 1 or 2")
  if(length(k) == 2 & min(k) == 1) stop("The minimum value of clusters to look for have to be 2")

  step <- match.arg(step)
  cut.method <- match.arg(cut.method)
  if(is.null(Beta.lower) & is.null(Beta.lower) & cut.method == "conf.int") stop("cut.method `conf.int' can not be used without Beta.lower and Beta.upper")
  if(cut.method == "conf.int"){
    code <- 0
    if((nrow(Beta.lower) != n) | (nrow(Beta.upper) != n)) code <- 1
    if((ncol(Beta.lower) != q) | (ncol(Beta.upper) != q)) code <- 2
    if(code != 0) stop("Dimensions of the matrices mismatched!")
  }
  if(alpha < 0 | alpha > 1) stop("alpha must be in (0,1)!")

  METHODS <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
  i.meth <- pmatch(method, METHODS)
  if(is.na(i.meth)) stop("invalid clustering method", paste("", method), "\n use one of this methods:", paste("", METHODS))
  if(i.meth == -1) stop("ambiguous clustering method", paste("", method))

  mat <- if(missing(diss.mat)) distshape(Beta, alpha = alpha, step = step, trace = trace) else diss.mat
  ogg <- hclust(mat, method=method)

  if(plot){
    par(mar = c(4,3.3,1.5,1.5)+.1, mfrow = c(1, 3))
    plot(ogg, main="", xlab="", ylab="", ann=TRUE, sub="")
    mtext(side=2, line=2, text="Height")
  }

  distance <- double(max(1, length(k)))

  if(length(k) == 1){
    pick <- k
  }
  else{
    if(ask){
      par(mfrow=c(1, 1))
      plot(ogg, main="", xlab="", ylab="", ann=TRUE, sub="")
      mtext(side=2, line=2, text="Height")

      pick <- 1
      while(pick > 0 && pick < (q+1)){
        pick <- menu(1:q, title = "\n Select how many clusters (or 0 to exit):\n")
        if(pick > 0 && pick < (q+1)){
          break
        }
        if(pick == 0) stop("Cluster algorithm aborted!!!")
      }
    }
    else{
      if(trace) cat("\nChoosing the optimal number of clusters: ", sep="")
      pick <- switch(cut.method,
                     "mindist"={
                       for(pick in k[1]:k[2]){
                         oc <- cutree(ogg, k = pick)
                         distance[pick - 1] <- .get_ave_sil_width(mat, oc)
                       }
                       names(distance) <- k[1]:k[2]
                       (k[1]:k[2])[which.max(distance)]},
                     "length"={
                       distance <- diff(rev(ogg$height))
                       names(distance) <- 2:(length(distance)+1)
                       tempIndex <- which(distance < 0)
                       tempIndex2 <- which.min(distance[tempIndex])
                       as.integer(names(tempIndex2))},
                     "conf.int"={
                       for(pick in k[1]:k[2]){
                         oc <- cutree(ogg, k = pick)
                         num.clust <- length(unique(oc))
                         tempDist <- vector(length = num.clust)
                         for(.i in unique(oc)){
                           indCl <- (oc == .i)
                           BetaLmean <- rowMeans(Beta.lower[, indCl, drop = FALSE])
                           BetaRmean <- rowMeans(Beta.upper[, indCl, drop = FALSE])

                           xx1 <- Beta[, indCl, drop = FALSE] - BetaRmean
                           xx2 <- Beta[, indCl, drop = FALSE] - BetaLmean

                           tempDist[.i] <- mean((colMeans(xx1 < 0) + colMeans(xx2 >= 0)) / 2)
                         }
                         distance[pick - 1] <- sum(tempDist * summary(silhouette(oc, mat))$clus.avg.widths) / sum(tempDist)
                       }
                       names(distance) <- k[1]:k[2]
                       (k[1]:k[2])[which.max(distance)]})
      if(trace) cat(paste0(pick, "\n\n"), sep="")
    }
  }

  if(plot) abline(h = mean(c(max(ogg$height), rev(ogg$height), 0)[c(pick, (pick+1))]), col=2, lty=2)

  oc <- cutree(ogg, k = pick)
  num.clust <- length(unique(oc))

  tab <- table(oc); ord <- order(tab); oc2 <- rep(0, length(oc))
  for(i in 1:num.clust) oc2 <- replace(oc2, oc == ord[i], values = i)
  oc <- oc2

  BetaUpper <- BetaLower <- BetaMedio <- Signif.interval <- NULL
  for(i in 1:num.clust){
    BetaMedio <- cbind(BetaMedio, rowMeans(Beta[, oc == i, drop = FALSE]))
    if(!is.null(Beta.lower) & !is.null(Beta.upper)){
      BetaLower <- cbind(BetaLower, rowMeans(Beta.lower[, oc == i, drop = FALSE]))
      BetaUpper <- cbind(BetaUpper, rowMeans(Beta.upper[, oc == i, drop = FALSE]))
      Signif.interval <- cbind(Signif.interval, apply(cbind(BetaLower[, i], BetaUpper[, i]), 1, prod))
    }
    else{
      a <- (1 - conf.level)/2; a <- c(a, 1 - a)
      BetaLower <- cbind(BetaLower, apply(Beta[, oc == i, drop = FALSE], 1, quantile, probs = .1))
      BetaUpper <- cbind(BetaUpper, apply(Beta[, oc == i, drop = FALSE], 1, quantile, probs = .9))
    }
    # if(plot){
    #   if(i == 1) {
    #     matplot(p, Beta[, oc == i, drop = FALSE], type="l", col=i, lty=2, ylim=rangeBeta, lwd=.8, xlab="", ylab="", axes=FALSE)
    #     axis(1); axis(2); mtext(side=1, line=2, text="p"); mtext(side=2, line=2, text="s(p)")
    #     lines(p, BetaMedio[, 1], col=1, lwd=1)
    #   }
    #   else{
    #     matlines(p, Beta[, which(oc == i), drop = FALSE], col=i, lwd=.8, lty=2)
    #     lines(p, BetaMedio[, i], col=i, lwd=1)
    #   }
    #   if(!is.null(Beta.lower) & !is.null(Beta.upper)){
    #     yy <- c(BetaLower[, i], tail(BetaUpper[, i], 1), rev(BetaUpper[, i]), BetaLower[1, i])
    #     xx <- c(p, tail(p, 1), rev(p), p[1])
    #     polygon(xx, yy, col = adjustcolor(i, alpha.f = 0.25), border = NA)
    #   }
    # }
  }

  # BetaDist <- sapply(1:num.clust, function(.i) (Beta[, oc == .i, drop = FALSE] - BetaMedio[, .i])^2, simplify=FALSE)
  # BetaDistMedio <- lapply(BetaDist, function(.x) sqrt(colSums(.x)))

  BetaDistMedio <- sapply(1:num.clust, function(.i) c(sqrt(2 * (1 - suppressWarnings(cor(Beta[, oc == .i, drop = FALSE], BetaMedio[, .i]))))), simplify = FALSE)

  oggSil <- if(length(unique(oc)) > 1) silhouette(oc, mat) else NULL

  DissMedio <- sapply(seq_len(num.clust), function(.i){
    tt <- which(oc == .i)
    tempcont <- NULL
    if(length(tt) == 1){
      tempcont <- c(tempcont, 0)
    }
    else{
      tempcont <- as.matrix(mat)[tt,tt][upper.tri(as.matrix(mat)[tt,tt])]
    }
    tempcont
  }, simplify=FALSE)

  names(BetaDistMedio) <- names(DissMedio) <- seq_len(num.clust) #names(BetaDist) <-

  if(plot){
    if(length(BetaDistMedio) > 0){
      boxplot(split(oggSil[,3], oc), names=as.integer(names(table(oc))), ylab="", axes=FALSE)
      axis(1, at=1:num.clust, labels=1:num.clust); axis(2); mtext(side=2, line=2, text="Silhouette")
    }
    if(ncol(BetaMedio) > 1){
      rangeBeta <- range(cbind(BetaMedio, BetaLower, BetaUpper))
      matplot(p, BetaMedio, type="l", lty=1, lwd=1, axes=FALSE, xlab="", ylab="", ylim=rangeBeta)
      axis(1); axis(2); mtext(side=1, line=2, text="p"); mtext(side=2, line=2, text="s(p)")
      for(i in 1:num.clust){
        yy <- c(BetaLower[, i], tail(BetaUpper[, i], 1), rev(BetaUpper[, i]), BetaLower[1, i])
        xx <- c(p, tail(p, 1), rev(p), p[1])
        polygon(xx, yy, col = adjustcolor(i, alpha.f = 0.25), border = NA)
      }
    }
    par(mfrow=c(1,1), mar=c(5,4,4,2)+.1)
  }

  obj <- list(call = match.call(), p = p, X = Beta, clusters = oc,
              X.mean = BetaMedio, X.mean.dist = BetaDistMedio,
              X.lower = Beta.lower, X.mean.lower = BetaLower,
              X.upper = Beta.upper, X.mean.upper = BetaUpper,
              Signif.interval = (Signif.interval > 0), k = pick,
              diss.matrix = mat, X.mean.diss = DissMedio,
              oggSilhouette = oggSil, oggHclust = ogg,
              distance = distance, step = step, method = method, cut.method = cut.method, alpha = alpha)
  class(obj) <- "clustEff"

  return(obj)
}

#' @export
distshape <- function(Beta, alpha=.5, step=c("both", "shape", "distance"), trace=TRUE){
  if(!is.matrix(Beta)) Beta <- as.matrix(Beta)
  n <- nrow(Beta)
  q <- ncol(Beta)
  p <- 1:n/n
  step <- match.arg(step)
  id <- combn(1:q, 2)
  mat <- matrix(NA, q, q)
  matShape <- matDist <- matrix(NA, n, ncol(id), dimnames = list(1:n, paste0(id[1,], "-", id[2,])))
  alphaDist <- NULL
  tem <- tem2 <- TRUE

  if(trace) cat("\nComputing Dissimilarity matrix: \n", sep="")
  if(step %in% c("both", "shape")){
    Smat <- apply(Beta, 2, splinefun, x=p)
    Smat2 <- sign(sapply(1:length(Smat), function(.i) Smat[[.i]](1:n/n, deriv=1)))
  }
  if(trace) pb <- txtProgressBar(min=0, max=ncol(id), style=3)
  for(i in 1:ncol(id)){
    if(trace) setTxtProgressBar(pb, i)
    if(step %in% c("both", "distance")){
      matDist[, i] <- abs(Beta[, id[1, i]] - Beta[, id[2, i]])
      # matDist[i, ] <- sqrt(2*(1-cor(s11, s22)))
    }
    if(step %in% c("both", "shape")){
      matShape[,i] <- (Smat2[, id[1, i]] * Smat2[, id[2, i]] == 1)
    }
  }
  if(step %in% c("both", "distance")){
    alphaDist <- apply(matDist, 1, function(.x) quantile(.x, alpha))
    matDist <- (matDist <= alphaDist)
  }
  if(trace) close(pb)

  mat[t(id)] <- .5*(colMeans(matDist) + colMeans(matShape))
  mat <- 1 - as.dist(t(mat))

  return(mat)
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

  avClust <- unlist(lapply(object$X.mean.dist, mean, na.rm=TRUE))
  avDiss <- unlist(lapply(object$X.mean.diss, mean))
  avSilhouette <- summary(object$oggSilhouette)$clus.avg.widths
  if(any(nClust == 1)) avSilhouette[nClust == 1] <- 1
  avClust <- c(avClust, sum(avClust * prop.table(nClust)))
  avDiss <- c(avDiss, sum(avDiss * prop.table(nClust)))
  avSilhouette <- c(avSilhouette, sum(avSilhouette * prop.table(nClust)))
  names(avSilhouette) <- names(avDiss) <- names(avClust) <- c(paste0("Cluster", seq_len(object$k)), "Average")

  obj <- list(call = object$call, k = object$k, n = nrow(object$X), p = ncol(object$X),
              step = object$step, alpha = object$alpha, method = object$method,
              cut.method = object$cut.method, tabClust = table(object$clusters),
              avClust = avClust, avSilhouette = avSilhouette, avDiss = avDiss)

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
plot.clustEff <- function(x, xvar=c("clusters", "dendrogram", "boxplot", "numclust"), which,
                          polygon=TRUE, dissimilarity=TRUE, par=FALSE, ...){
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
           if(is.null(L$ylab)) L$ylab <- "Height"
           if(is.null(L$lty)) L$lty <- 2
           if(is.null(L$mar)) L$mar <- c(2,3,4,1)+.1
           par(mar=L$mar)

           plot(x$oggHclust, main=L$main, xlab="", ylab="", ann=TRUE, sub="")
           mtext(side=1, line=2, text=L$xlab)
           mtext(side=2, line=2, text=L$ylab)
           abline(h=mean(c(max(x$oggHclust$height), rev(x$oggHclust$height), 0)[c(x$k, (x$k+1))]), col=L$col, lty=L$lty)
         },
         "clusters"={
           if(missing(which)) which <- seq(x$k)
           yRange <- range(x$X.mean[,  which])
           if(!is.null(x$X.mean.lower) & !is.null(x$X.mean.upper)){
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
           if(is.null(L$lty)){L$lty <- c(1, 3)}
           if(is.null(L$col)) L$col <- gray.colors(length(which))
           if(length(L$col) == 1) L$col <- rep(L$col, length(which))

           if(length(which) == 1){
             if(is.null(L$mar)) L$mar <- c(4,3.3,2,1)+.1
             par(mar=L$mar)

             tempInd <- x$clusters == which
             matplot(x$p, x$X[, tempInd], xlab="", ylab="", type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty[2], axes=FALSE)
             axis(1); axis(2); mtext(side=1, line=2, text=L$xlab); mtext(side=2, line=2, text=L$ylab)
             lines(x$p, x$X.mean[, which], lwd=L$lwd[2], col=L$col[1], lty = L$lty[1])

             if(!is.null(x$X.mean.lower) & !is.null(x$X.mean.upper)){
               if(!is.null(polygon)){
                 if(polygon){
                   yy <- c(x$X.mean.lower[, which], tail(x$X.mean.upper[, which], 1), rev(x$X.mean.upper[, which]), x$X.mean.lower[1, which])
                   xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
                   polygon(xx, yy, col = adjustcolor(L$col[1], alpha.f = .25), border = NA)
                 }
                 else{
                   points(x$p, x$X.mean.lower[, which], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[1])
                   points(x$p, x$X.mean.upper[, which], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[1])
                 }
                 matlines(x$p, x$X[, tempInd], xlab="", ylab="", type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty[2], axes=FALSE)
                 lines(x$p, x$X.mean[, which], lwd=L$lwd[2], col=L$col[1], lty = L$lty[1])
               }
             }
             abline(h=0, lty=3, col=1)
           }
           else{
             if(is.null(L$mar)) L$mar <- c(4,3.3,1,1)+.1
             par(mar=L$mar)

             matplot(x$p, x$X.mean, xlab="", ylab="", type=L$type, ylim=L$ylim,
                     lwd=L$lwd[1], main=L$main, col=L$col, lty=L$lty[1], axes=FALSE)
             # matplot(x$p, x$X[, x$clusters == which[1]], xlab="", ylab="", type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty, axes=FALSE)
             axis(1); axis(2); mtext(side=1, line=2, text=L$xlab); mtext(side=2, line=2, text=L$ylab)

             if(!is.null(x$X.mean.lower) & !is.null(x$X.mean.upper)){
               if(!is.null(polygon)){
                 for(i in which){
                   if(polygon){
                     yy <- c(x$X.mean.lower[, i], tail(x$X.mean.upper[, i], 1), rev(x$X.mean.upper[, i]), x$X.mean.lower[1, i])
                     xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
                     polygon(xx, yy, col = adjustcolor(L$col[i], alpha.f = .25), border = NA)# adjustcolor(L$col[1], alpha.f = 0.25)
                   }else{
                     points(x$p, x$X.mean.lower[, i], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[i])
                     points(x$p, x$X.mean.upper[, i], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[i])
                   }
                 }
                 matlines(x$p, x$X.mean, xlab="", ylab="", type=L$type, ylim=L$ylim,
                          lwd=L$lwd[2], main=L$main, col=L$col, lty=L$lty[1], axes=FALSE)
               }
             }
           }
         },
         "boxplot"={
           tabClust <- table(x$clusters)
           k <- length(tabClust)
           X <- if(dissimilarity) x$X.mean.diss else split(x$oggSilhouette[, 3], x$clusters)

           if(is.null(L$main)) L$main <-  if(dissimilarity) "Average dissimilarity within cluster" else "Silhouette within cluster"
           if(is.null(L$labels)) L$labels <- as.integer(names(tabClust))[tabClust > 1]
           if(is.null(L$ylab)) L$ylab <- if(dissimilarity) "Dissimilarity" else "Silhouette"
           if(is.null(L$ylim)) L$ylim <- c(0, 1)
           if(!dissimilarity) L$ylim <- c(0, max(unlist(X)))
           if(is.null(L$mar)) L$mar <- c(3.3,3.3,3.3,1)+.1
           par(mar=L$mar)

           # par(mfrow=c(1,1))
           if(length(X) > 0){
             boxplot(X, names=L$labels, ylim=L$ylim, main=L$main, ylab="", axes=FALSE)
             axis(1, at=1:k, labels=1:k); axis(2); mtext(side=2, line=2, text=L$ylab)
           }
         },
         "numclust"={
           if(is.null(L$xlab)) L$xlab <- "N. of clusters"
           if(is.null(L$ylab)) L$ylab <- expression(pi[out]^k-pi[out]^(k+1))
           if(is.null(L$main)) L$main <- ""
           if(is.null(L$lwd)) L$lwd <- 1
           if(is.null(L$pch)) L$pch <- 19
           if(is.null(L$type)) L$type <- "b"
           if(is.null(L$col)) L$col <- 1
           if(is.null(L$mar)) L$mar <- c(4,4,2,1)+.1
           par(mar=L$mar)

           plot(diff(x$distance), type=L$type, axes=FALSE, col=L$col, lwd=L$lwd, xlab="", ylab="", pch=L$pch)
           axis(1, at=1:length(diff(x$distance)), labels=names(x$distance)[-1])
           axis(2)
           mtext(side=1, line=2, text=L$xlab)
           mtext(side=2, line=2, text=L$ylab)
         })

  par(mar=c(5,4,4,2)+.1)
}

#' @export
extract.object <- function(Y, X, intercept=TRUE, formula.p=~slp(p, 3), s, object, p, which){
  if(missing(p)) p <- seq(.01, .99, .01)

  if(!missing(object)){
    if(!inherits(object, "iqr")){
      stop("Wrong class object!")
    }
    else{
      if(is.data.frame(object$mf)){
        X <- as.matrix(object$mf[, -1])
      }else{
        X <- object$mf[[2]]
      }

      n <- nrow(X)
      q <- nrow(object$coefficients)
      intercept <- attr(attr(object$mf, "terms"), "intercept")
      if(missing(which)) which <- 1:q
      labels <- rownames(object$coefficients)

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
      labels <- paste0(X, 1:(qX - 1))
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
    }
    else{
      B <- model.matrix(formula.p, data = data.frame(p=p))
      k <- ncol(B)
      termlabelsB <- attr(terms(formula.p), "term.labels")
      coefnamesB <- colnames(B)
    }
    if(missing(s)) s <- matrix(1, nrow=qX, ncol=k)
    colnames(s) <- coefnamesB
    rownames(s) <- labels
    if(missing(which)) which <- 1:(qX-intercept)

    if(qY == 1){
      object <- if(intercept) iqr(Y[, 1] ~ X, formula.p=formula.p, s=s) else iqr(Y[, 1] ~ -1 + X, formula.p=formula.p, s=s)
      predObj <- predict(object, type="beta", p=p)
      tempX <- sapply((2-!intercept):qX, function(i) predObj[[i]][,2])[, which]
      tempXl <- sapply((2-!intercept):qX, function(i) predObj[[i]][,4])[, which]
      tempXr <- sapply((2-!intercept):qX, function(i) predObj[[i]][,5])[, which]
      colnames(tempX) <- colnames(tempXl) <- colnames(tempXr) <- if(intercept) labels[-1][which] else labels[which]
    }
    else{
      predObj <- list()
      for(i in 1:qY){
        object <- if(intercept) iqr(Y[, i] ~ X, formula.p=formula.p, s=s) else iqr(Y[, i] ~ -1 + X, formula.p=formula.p, s=s)
        predObj[[i]] <- predict(object, type="beta", p=p)
      }
      predObjX <- lapply(predObj, function(.x) simplify2array(lapply(.x, function(.x2) .x2[,2]))[, which+intercept, drop=FALSE])
      predObjXl <- lapply(predObj, function(.x) simplify2array(lapply(.x, function(.x2) .x2[,4]))[, which+intercept, drop=FALSE])
      predObjXr <- lapply(predObj, function(.x) simplify2array(lapply(.x, function(.x2) .x2[,5]))[, which+intercept, drop=FALSE])

      tempX <- tempXl <- tempXr <- list()
      for(i in 1:ncol(predObjX[[1]])){
        tempX[[i]] <- sapply(1:qY, function(j) predObjX[[j]][, i])
        tempXl[[i]] <- sapply(1:qY, function(j) predObjXl[[j]][, i])
        tempXr[[i]] <- sapply(1:qY, function(j) predObjXr[[j]][, i])
      }
      names(tempX) <- names(tempXl) <- names(tempXr) <- if(intercept) labels[-1][which] else labels[which]
    }
  }

  return(list(p=p, X=tempX, Xl=tempXl, Xr=tempXr))
}


# extract.object.interaction <- function(obj, extr_obj, p=1:99/100, interaction){
#   n_list <- length(interaction)
#   n_obj <- length(extr_obj$X)
#   temp_temp <- summary(obj, p=p, cov=T)
#   for(i in 1:n_list){
#     temp_int <- interaction[[i]]
#     se1 <- sapply(1:length(p), function(.i){
#       temp_cov <- temp_temp[[.i]]$cov[temp_int, temp_int]
#       sum(diag(temp_cov)) - 2*temp_cov[1,2]
#     })
#     extr_obj$X$interaction <- extr_obj$X[[temp_int[1]]]+ extr_obj$X[[temp_int[2]]]
#     extr_obj$Xl$interaction <- extr_obj$X$interaction - 1.96*se1
#     extr_obj$Xr$interaction <- extr_obj$X$interaction + 1.96*se1
#     names(extr_obj$X)[n_obj+i] <- names(extr_obj$Xl)[n_obj+i] <- names(extr_obj$Xr)[n_obj+i] <- paste0("interaction",i)
#   }
#   return(extr_obj)
# }

#' @export
fpcac <- function(X, K = 2, fd = NULL, nbasis = 5, norder = 3, nharmonics = 3,
                  alpha = 0, niter = 30, Ksteps = 25, conf.level = 0.9, seed, disp = FALSE){

  if(is.null(fd)) {
    Xorig <- as.matrix(X)
    bspl <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
    fd <- Data2fd(argvals = seq(0, 1, length.out = nrow(Xorig)), y = Xorig, basisobj = fd(basisobj = bspl))
  }
  else{
    Xorig <- eval.fd(fd[[3]][[1]]/max(fd[[3]][[1]]), fd, Lfdobj=0)
  }
  X <- pca.fd(fd, nharm = nharmonics)$scores

  n <- dim(X)[1]
  p <- dim(X)[2]
  no.trim <- floor(n*(1-alpha))
  ll <- double(K)
  ind <- integer(n)
  dist <- ind
  seqK <- seq.int(K)
  seqN <- seq.int(n)
  seqNoTrim <- seq.int(no.trim)

  # nit = Number or random restarting (larger values provide more accurate solutions
  # Ksteps = Number of k-Mean steps (not too many ksteps are needed)

  # Initialize the objective function by a large enough value
  vopt <- 1e+6

  # set random seed if not missing
  if(!missing(seed)) set.seed(seed)

  #Ramdon restarts
  for (iter in 1:niter){
    # Randomly choose the K initial centers
    cini <- X[sample(n, size=K, replace=F), ]
    dim(cini) <- c(K, p)

    # C-steps
    for (t in 1:Ksteps){
      # Distances of each data point to its closest center
      ll <- sapply(seqK, function(k) rowSums((X - rep(cini[k,], rep.int(n, p)))^2))
      dist <- apply(ll, 1, min)
      ind <- apply(ll, 1, which.min)

      # Modified data (Xmod) with the non-trimmed points and last
      # column equal to the clusters allocations
      qq <- seqN[dist <= sort(dist)[no.trim]]
      xmod <- cbind(X[qq, ], ind[qq])

      #Calculus of the new k centers
      cini <- t(sapply(seqK, function(k) apply(xmod[xmod[, p+1] == k, 1:p, drop=FALSE], 2, mean, na.rm=TRUE)))
    }
    if(any(is.na(cini))) next

    # Calculus of the trimmed k-variance
    obj <- mean(sapply(seqNoTrim, function(l) sum((xmod[l, 1:p] - cini[xmod[l, p+1], ])^2)))

    # Change the optimal value and the optimal centers (copt)
    # if a reduction in the objective function happens
    if (obj < vopt){
      vopt <- obj
      # Improvements in the objective functions are printed
      if(disp){
        cat("\n")
        cat("Iter ", iter, "\t Objective function = ",
            formatC(vopt, digits=5, width=8, format="f"),
            "<----")
      }
      copt <- cini
    }
    else{
      if(disp){
        cat("\n")
        cat("Iter ", iter, "\t Objective function = ", formatC(obj, digits=5, width=8, format="f"))
      }
    }
  }

  ## Obtain the final cluster allocations (this is necesary, because a final
  # cluster assigment movement is possible)
  asig <- ind
  ll <- sapply(seqK, function(k) rowSums((X - rep(copt[k,], rep.int(n, p)))^2))
  dist <- apply(ll, 1, min)
  ind <- apply(ll, 1, which.min)

  # Compute the radius of the optimal ball ropt
  ord <- sort(dist)
  ropt <- ord[no.trim]

  # Assign every observation to each cluster and 0 for the trimmed observations
  asig <- ifelse(dist > ropt, 0, ind)

  tab <- table(asig); if(length(tab) > K) {tab <- tab[-1]}
  ord <- order(tab); copt <- copt[ord, ]; asig2 <- rep(0, length(asig))
  for(i in 1:K) asig2 <- replace(asig2, asig == ord[i], values = i)
  asig <- asig2

  # Print the clusters results
  # if(disp){
  #   cat("\n")
  #   for (k in seqK){
  #     # Group
  #     cat("\nCluster ", k, ":\n")
  #     print(seqN[asig == k])
  #   }
  #   # Trimmed observations
  #   if(sum(asig == 0) == 0)
  #     cat("\nNo trimmed observations!")
  #   else{
  #     cat("\nTrimmed observations:")
  #     print(seqN[asig == 0])
  #   }
  # }

  # calculate mean curves and silhouette
  X.mean <- sapply(1:K, function(i) rowMeans(Xorig[, asig == i, drop = FALSE]))
  X.mean.dist <- sapply(1:K, function(.i) c(sqrt(2 * (1 - suppressWarnings(cor(Xorig[, asig == .i, drop = FALSE], X.mean[, .i]))))), simplify = FALSE)
  a <- (1 - conf.level)/2; a <- c(a, 1 - a)
  X.mean.lower <- sapply(1:K, function(i) apply(Xorig[, asig == i, drop = FALSE], 1, quantile, probs = a[1]))
  X.mean.upper <- sapply(1:K, function(i) apply(Xorig[, asig == i, drop = FALSE], 1, quantile, probs = a[2]))
  # diss <- dist(X)
  # oggSil <- silhouette(asig, diss)

  # Function trimm returns the objetive value (vopt), the trimmed k-means centers (copt),
  # the optimal radius (ropt) and the assignation of each observation (asig=0 are the
  # trimmed ones
  ris <- list("obj.function" = vopt, "centers" = copt, "radius" = ropt, "clusters" = asig,
              "Xorig" = Xorig, "fd" = fd, "X" = X , "X.mean" = X.mean, "X.mean.dist" = X.mean.dist,
              "X.mean.lower" = X.mean.lower, "X.mean.upper" = X.mean.upper,
              # "diss.matrix" = diss, "oggSilhouette" = oggSil,
              "call"=match.call())
  class(ris) <- "fpcac"

  return(ris)
}

#' @export
print.fpcac <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  tabClust <- as.integer(table(x$clusters[x$clusters > 0]))
  trimmed <- sum(x$clusters == 0)

  cat("###################################", "\n")
  cat("n. of clusters:", format(nrow(x$centers), digits=digits), "\n")
  cat("n. of curves:", length(x$clusters), "\n")
  cat("n. of harmonics:", ncol(x$centers), "\n")
  cat("n. of trimmed curves:", trimmed, "\n")
  cat("###################################", "\n\n")

  cat("FPCAC algorithm with ", nrow(x$centers), " clusters of sizes ",
      paste(tabClust, collapse = ", "), "\n", sep = "")
  cat("\n")

  invisible(x)
}

#' @export
summary.fpcac <- function(object, ...){
  tabClust <- table(object$clusters)
  trimmed <- sum(object$clusters == 0)
  if(trimmed > 0) tabClust <- tabClust[-1]
  nClust <- length(tabClust)

  avClust <- unlist(lapply(object$X.mean.dist, mean, na.rm=TRUE))
  # if(trimmed > 0) avClust <- avClust[-1]
  # if(any(nClust == 1)) avClust[nClust == 1] <- 1
  avClust <- c(avClust, sum(avClust * prop.table(tabClust)))
  names(avClust) <- c(paste0("Cluster", seq_len(nClust)), "Average")

  obj <- list(call = object$call, k = nClust, n = nrow(object$X), p = ncol(object$X),
              trimmed = trimmed, tabClust = tabClust, avClust = avClust)

  class(obj) <- "summary.fpcac"
  return(obj)
}

#' @export
print.summary.fpcac <- function(x, digits=max(3L, getOption("digits") - 3L), ...){

  cat("\n######################", "\n\n")
  cat("Selected n. of clusters:", format(x$k, digits=digits), "\n")
  cat("n. of curves:", x$n, "\n")
  cat("n. of harmonics:", x$p, "\n")
  cat("n. of trimmed curves:", x$trimmed, "\n")
  cat("\n######################", "\n\n")

  cat("Clustering table:")
  print(x$tabClust, ...)
  cat("\n######################", "\n\n")

  cat("Average cluster distance:\n")
  print(round(x$avClust, digits), ...)
  cat("\n######################", "\n\n")

  # if(!is.null(x$avSilhouette)){
  #   cat("Individual silhouette widths:\n")
  #   print(round(x$avSilhouette, digits), ...)
  #   cat("\n######################", "\n\n")
  # }

  invisible(x)
}

#' @export
plot.fpcac <- function(x, which, polygon=TRUE, conf.level, ...){
  if(!missing(which)){
    if(any(which <= 0) | any(which > x$k)){
      stop("Which values in 1-", x$k)
    }
  }

  x$k <- nrow(x$centers)
  x$p <- seq_len(nrow(x$X.mean))
  L <- list(...)

  if(!missing(conf.level)){
    a <- (1 - conf.level)/2; a <- c(a, 1 - a)
    x$X.mean.lower <- sapply(1:x$k, function(i) apply(x$Xorig[, x$clusters == i, drop = FALSE], 1, quantile, probs = a[1]))
    x$X.mean.upper <- sapply(1:x$k, function(i) apply(x$Xorig[, x$clusters == i, drop = FALSE], 1, quantile, probs = a[2]))
  }

  if(missing(which)) which <- seq(x$k)
  yRange <- range(x$X.mean[,  which])
  temp <- range(cbind(x$X.mean.lower, x$X.mean.upper))
  yRange <- range(c(temp, yRange))

  tabClust <- table(x$clusters)

  if(is.null(L$xlab)) L$xlab <- "Time"
  if(is.null(L$ylab)) L$ylab <- "F(time)"
  if(is.null(L$main)) L$main <- ""
  if(is.null(L$lwd)) L$lwd <- c(1, 1.5)
  if(length(L$lwd) == 1) L$lwd <- rep(L$lwd, 2)
  if(is.null(L$type)) L$type <- "l"
  if(is.null(L$ylim)) L$ylim <- yRange
  if(is.null(L$lty)){L$lty <- c(1, 3)}
  if(is.null(L$col)) L$col <- gray.colors(length(which))
  if(length(L$col) == 1) L$col <- rep(L$col, length(which))

  if(length(which) == 1){
    if(is.null(L$mar)) L$mar <- c(4,3.3,2,1)+.1
    par(mar=L$mar)

    tempInd <- x$clusters == which
    matplot(x$p, x$Xorig[, tempInd, drop = FALSE], xlab="", ylab="", type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty[2], axes=FALSE)
    axis(1); axis(2); mtext(side=1, line=2, text=L$xlab); mtext(side=2, line=2, text=L$ylab)
    lines(x$p, x$X.mean[, which], lwd=L$lwd[2], col=L$col[1], lty = L$lty[1])

    if(!is.null(polygon)){
      if(polygon){
        yy <- c(x$X.mean.lower[, which], tail(x$X.mean.upper[, which], 1), rev(x$X.mean.upper[, which]), x$X.mean.lower[1, which])
        xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
        polygon(xx, yy, col = adjustcolor(L$col[1], alpha.f = .25), border = NA)
      }
      else{
        points(x$p, x$X.mean.lower[, which], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[1])
        points(x$p, x$X.mean.upper[, which], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[1])
      }
      matlines(x$p, x$Xorig[, tempInd, drop = FALSE], xlab="", ylab="", type=L$type, ylim=L$ylim, lwd=L$lwd[1], main=L$main, col=L$col[1], lty=L$lty[2], axes=FALSE)
      lines(x$p, x$X.mean[, which], lwd=L$lwd[2], col=L$col[1], lty = L$lty[1])
    }

    abline(h=0, lty=3, col=1)
  }
  else{
    if(is.null(L$mar)) L$mar <- c(4,3.3,1,1)+.1
    par(mar=L$mar)

    matplot(x$p, x$X.mean, xlab="", ylab="", type=L$type, ylim=L$ylim,
            lwd=L$lwd[1], main=L$main, col=L$col, lty=L$lty[1], axes=FALSE)
    axis(1); axis(2); mtext(side=1, line=2, text=L$xlab); mtext(side=2, line=2, text=L$ylab)

    if(!is.null(polygon)){
      for(i in which){
        if(polygon){
          yy <- c(x$X.mean.lower[, i], tail(x$X.mean.upper[, i], 1), rev(x$X.mean.upper[, i]), x$X.mean.lower[1, i])
          xx <- c(x$p, tail(x$p, 1), rev(x$p), x$p[1])
          polygon(xx, yy, col = adjustcolor(L$col[i], alpha.f = .25), border = NA)# adjustcolor(L$col[1], alpha.f = 0.25)
        }else{
          points(x$p, x$X.mean.lower[, i], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[i])
          points(x$p, x$X.mean.upper[, i], lty = 2, lwd = L$lwd[2], type = "l", col = L$col[i])
        }
      }
      matlines(x$p, x$X.mean, xlab="", ylab="", type=L$type, ylim=L$ylim,
               lwd=L$lwd[2], main=L$main, col=L$col, lty=L$lty[1], axes=FALSE)
    }
  }

  par(mar=c(5,4,4,2)+.1)
}

#' @export
opt.fpcac <- function(X, k.max = 5, method = c("silhouette", "wss"),
                      fd = NULL, nbasis = 5, norder = 3, nharmonics = 3,
                      alpha = 0, niter = 30, Ksteps = 10, seed,
                      diss = NULL, trace=FALSE){
  if(!missing(seed)) set.seed(seed)
  method <- match.arg(method)

  if(trace) pb <- txtProgressBar(min=0, max=k.max, style=3)
  objTemp <- sapply(seq_len(k.max), function(.k) {
    if(trace) setTxtProgressBar(pb, .k)
    if(is.null(fd)) fpcac(X, K=.k, alpha=alpha, niter=niter, Ksteps=Ksteps, nbasis=nbasis, norder=norder, nharmonics=nharmonics, disp=FALSE)
    else fpcac(fd=fd, K=.k, alpha=alpha, niter=niter, Ksteps=Ksteps, nbasis=nbasis, norder=norder, nharmonics=nharmonics, disp=FALSE)
  }, simplify = FALSE)
  if(trace) close(pb)
  names(objTemp) <- seq_len(k.max)

  clusters <- simplify2array(lapply(objTemp, function(.x) .x$clusters))

  nfolds <- dim(objTemp[[1]]$X)[1] %/% 5000
  nfolds <- pmax(nfolds, 1)
  foldid <- sort(rep(seq(nfolds), length = dim(objTemp[[1]]$X)[1]))
  if(is.null(diss)) {
    v <- NULL
    for(i in 1:nfolds) {
      # diss <- dist(objTemp[[1]]$X[foldid == i, ])
      diss <- as.dist(sqrt(2*(1-cor(t(objTemp[[1]]$X[foldid == i, ])))))
      temp <- if (method == "silhouette")
        c(0, apply(clusters[foldid == i, -1], 2, .get_ave_sil_width, d=diss))
      else
        apply(clusters, 2, .get_withinSS, d=diss)
      v <- rbind(v, temp)
    }
    v <- colMeans(v)
  }
  else{
    v <- if (method == "silhouette")
      c(0, apply(clusters[,-1], 2, .get_ave_sil_width, d=diss))
    else
      apply(clusters, 2, .get_withinSS, d=diss)
  }
  K.opt <- if(method == "wss") as.integer(names(which.min(diff(v)))) else which.max(v)

  df <- data.frame(clusters = as.factor(1:k.max), y = v)
  linecolor <- "steelblue"
  ylab <- "Total Within Sum of Square"
  if(method == "silhouette") ylab <- "Average silhouette width"
  p <- ggline(df, x = "clusters", y = "y", group = 1,
              color = linecolor, ylab = ylab, xlab = "Number of clusters k",
              main = "Optimal number of clusters")
  if(method == "silhouette")
    p <- p + geom_vline(xintercept = which.max(v), linetype = 2, color = linecolor)
  print(p)

  out <- list(obj.function=v, clusters=clusters, K=seq_len(k.max), K.opt=K.opt, plot=p)
  return(out)
}

.get_ave_sil_width <- function (d, cluster){
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("cluster package needed for this function to work. Please install it.")
  }
  ss <- silhouette(cluster, d)
  mean(ss[, 3])
}
.get_withinSS <- function (d, cluster){
  d <- as.dist(d)
  cn <- max(cluster)
  clusterf <- as.factor(cluster)
  clusterl <- levels(clusterf)
  cnn <- length(clusterl)
  if (cn != cnn) {
    warning("cluster renumbered because maximum != number of clusters")
    for (i in 1:cnn) cluster[clusterf == clusterl[i]] <- i
    cn <- cnn
  }
  cwn <- cn
  dmat <- as.matrix(d)
  within.cluster.ss <- 0
  for (i in 1:cn) {
    cluster.size <- sum(cluster == i)
    di <- as.dist(dmat[cluster == i, cluster == i])
    within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size
  }
  within.cluster.ss
}
