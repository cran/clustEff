\name{distshape}
\alias{distshape}
\title{
Dissimilarity matrix
}
\description{
This function implements the dissimilarity matrix based on shape and distance of curves.
}
\usage{
distshape(Beta, alpha=.5, step=c("both", "shape", "distance"), trace=TRUE)
}
\arguments{
  \item{Beta}{
    A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of curves to cluster and \kbd{n} is either the length of percentiles used in the quantile regression or the length of the time vector.
  }
  \item{alpha}{
    It is the alpha-percentile used for computing the dissimilarity matrix. If not fixed, the algorithm choose alpha=.25 (cluster.effects=TRUE) or alpha=.5 (cluster.effects=FALSE).
  }
  \item{step}{
    The steps used in computing the dissimilarity matrix. Default is "both"=("shape" and  "distance")
  }
  \item{trace}{
    If TRUE, some informations are printed.
  }
}
\value{
The dissimilarity matri of class \dQuote{\code{dist}}.
}
\references{
Sottile, G and Adelfio, G (2017). \emph{Clustering of effects through quantile regression}.
Proceedings 32nd International Workshop of Statistical Modeling, Groningen (NL), vol.2 127-130, https://iwsm2017.webhosting.rug.nl/IWSM_2017_V2.pdf.

Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, doi: 10.1111/biom.12410.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{clustEff}},\code{\link{summary.clustEff}}, \code{\link{plot.clustEff}},
for summary and plotting.
\code{\link{extract.object}} to extract useful objects for the clustering algorithm through a quantile regression coefficient modeling in a multivariate case.
}
\examples{

  ##### Using simulated data in all examples

  # see the documentation for 'clustEff-package'

}
\keyword{models}
\keyword{clustering algorithm}
