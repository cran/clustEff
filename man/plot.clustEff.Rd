\name{plot.clustEff}
\alias{plot.clustEff}
\title{
Plot Clustering Effects
}
\description{
Produces a dendrogram, a cluster plot and a boxplot of average distance cluster
class \dQuote{\code{piqr}}.
}
\usage{
\method{plot}{clustEff}(x, xvar=c("clusters", "dendrogram", "boxplot", "numclust"), which,
        add=FALSE, all=TRUE, polygon=TRUE, dissimilarity=TRUE, \ldots)
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
  \item{add}{If TRUE and xvar="clusters", a unique plot of clusters is done.}
  \item{all}{If TRUE and add=FALSE and xvar="clusters", a unique window of clusters is done.}
  \item{polygon}{If TRUE confidence intervals are represented by shaded areas via polygon. Otherwise, dashed lines are used.}
  \item{dissimilarity}{If TRUE dissimilarity measure within each cluster is used to do boxplot representation.}
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

  # see the documentation for 'clustEff-package'

}
\keyword{methods}

