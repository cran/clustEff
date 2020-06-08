\name{summary.fpcac}
\alias{summary.fpcac}
\title{
Summary FPCAC algorithm
}
\description{
Summary of an object of class \dQuote{\code{fpcac}}.
}
\usage{
\method{summary}{fpcac}(object, \ldots)
}
\arguments{
  \item{object}{
  An object of class \dQuote{\code{fpcac}}, the result of a call to \code{\link{fpcac}}.
}
  \item{\ldots}{for future methods.}
}
\details{
A summary of the clustering algorithm is printed.
}
\value{
The following items is returned:
  \item{k}{The number of selected clusters.}
  \item{n}{The number of curves.}
  \item{p}{The number of harmonics used.}
  \item{trimmed}{The number of trimmed curves.}
  \item{tabClust}{The table of clusters.}
  \item{avClust}{The average distance within clusters.}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{fpcac}}, \code{\link{opt.fpcac}}
}
\examples{

# using simulated data

# see the documentation for 'fpcac'

}
