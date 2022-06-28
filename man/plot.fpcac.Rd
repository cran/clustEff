\name{plot.fpcac}
\alias{plot.fpcac}
\title{
Plot Functional Principal Component Analysis Clustering
}
\description{
Produces a cluster plot of an object of class \dQuote{\code{fpcac}}.
}
\usage{
\method{plot}{fpcac}(x, which, polygon=TRUE, conf.level, \ldots)
}
\arguments{
  \item{x}{An object of class \dQuote{\code{clustEdd}}, typically the result of a call to \code{\link{fpcac}}.}
  \item{which}{If missing all curves effect are plotted.}
 \item{polygon}{If TRUE confidence intervals are represented by shaded areas via polygon. Otherwise, dashed lines are used. If NULL no confidence intervals are represented}
 \item{conf.level}{the confidence level required.}
  \item{\ldots}{additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd, lty}. See \code{\link{par}}.}
}
\details{
Different plot for the clustering algorithm.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.ot}
}
\seealso{
\code{\link{fpcac}}, \code{\link{summary.fpcac}}, \code{\link{opt.fpcac}}.
}
\examples{

  # using simulated data

  # see the documentation for 'fpcac'

}
\keyword{methods}

