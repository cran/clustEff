\name{summary.clustEff}
\alias{summary.clustEff}
\title{
Summary after Clustering Algorithm
}
\description{
Summary of an object of class \dQuote{\code{clustEff}}.
}
\usage{
\method{summary}{clustEff}(object, \ldots)
}
\arguments{
  \item{object}{
  An object of class \dQuote{\code{clustEff}}, the result of a call to \code{\link{clustEff}}.
}
  \item{\ldots}{for future methods.}
}
\details{
A summary of the clustering algorithm is printed.
}
\value{
The following items is returned:
  \item{k}{The number of selected clusters.}
  \item{n}{The number of observations.}
  \item{p}{The number of curves.}
  \item{step}{The selected step for computing the dissimilarity matrix.}
  \item{alpha}{The alpha-percentile used for computing the dissimilarity matrix.}
  \item{method}{The selected method to compute the hierarchical cluster analysis.}
  \item{cut.method}{The selected method to choose the best number of clusters.}
  \item{tabClust}{The table of clusters.}
  \item{avClust}{The average distance within clusters.}
  \item{avSilhouette}{Silhouette widths for clusters.}
  \item{avDiss}{The average dissimilarity measure within clusters.}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{clustEff}}, for cluster algorithm\code{\link{extract.object}} for extracting information through a quantile regression coefficient modeling in a multivariate case and plotting objects of class \dQuote{\code{clustEff}}.
}
\examples{

# using simulated data

# see the documentation for 'clustEff-package'

}
