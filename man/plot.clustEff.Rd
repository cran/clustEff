\name{plot.clustEff}
\alias{plot.clustEff}
\title{
Plot Clustering Effects
}
\description{
Produces a dendrogram, a cluster plot and a boxplot of average distance cluster of an object of class \dQuote{\code{clustEff}}.
}
\usage{
\method{plot}{clustEff}(x, xvar=c("clusters", "dendrogram", "boxplot", "numclust"), which,
        polygon=TRUE, dissimilarity=TRUE, par=FALSE, \ldots)
}
\arguments{
  \item{x}{An object of class \dQuote{\code{clustEdd}}, typically the result of a call to \code{\link{clustEff}}.}
\item{xvar}{
  Clusters: plot of the k clusters;
  Dendrogram: plot of the tree after computing the dissimilarity measure and applying a hierarchical       clustering algorithm;
  Boxplot: plot the average distance within clusters;
  Numclust: plot the curve to minimize to select the best number of clusters;
}
  \item{which}{If missing all curves effect are plotted.}
 \item{polygon}{If TRUE confidence intervals are represented by shaded areas via polygon. Otherwise, dashed lines are used. If NULL no confidence intervals are represented}
  \item{dissimilarity}{If TRUE dissimilarity measure within each cluster is used to do boxplot representation.}
  \item{par}{If TRUE the screen is automaticcaly splitted.}
  \item{\ldots}{additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd, lty}. See \code{\link{par}}.}
}
\details{
Different plot for the clustering algorithm.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.ot}
}
\seealso{
\code{\link{clustEff}} for cluster algorithm; \code{\link{extract.object}} for extracting information through a quantile regression coefficient modeling in a multivariate case; \code{\link{summary.clustEff}} for clustering summary.
}
\examples{

  # using simulated data

  # see the documentation for 'clustEff'

}
\keyword{methods}

