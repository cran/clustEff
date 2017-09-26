\name{extract.object}
\alias{extract.object}
\title{
\code{extract.object} fits a multivariate quantile regression and extracts objects for the cluster effects algorithm.
}
\description{
\code{extract.object} fits a multivariate quantile regression and extracts objects for the cluster effects algorithm.
}
\usage{
extract.object(Y, X, intercept=TRUE, formula.p=~slp(p, 3), s, object, p, which)
}
\arguments{
  \item{Y}{A multivariate response matrix of dimension \kbd{n} x \kbd{q1}, or a vector of length \kbd{n}.}
  \item{X}{The covariates matrix of dimension \kbd{n} x \kbd{q2}.}
  \item{intercept}{If TRUE, the intercept is included in the model.}
  \item{formula.p}{a one-sided formula of the form \code{~ b1(p, \ldots) + b2(p, \ldots) + \ldots}}
  \item{s}{An optional 0/1 matrix that allows to exclude some model coefficients (see \sQuote{Examples}).}
  \item{object}{An object of class \dQuote{\code{iqr}}. If missing, Y and X have to be supplied.}
  \item{p}{The percentiles used in quantile regression coefficient modeling. If missing a default sequence is choosen.}
  \item{which}{If fixed, only the selected covariates are extraced from the model. If missing all the covariates are extracted.}
}

\details{
A list of objects useful to run the cluster effect algorithm is created.
}

\value{
  \item{p}{The percentiles used in  the quantile regression.}
  \item{X}{A list containing as many matrices as covariates, where for each matrix the number of columns corresponds to the number of the responses. Each column of a matrix corresponds to one curve effect.}
  \item{Xl}{A list as \kbd{X}. Each column of a matrix corresponds to the lower interval of the curve effect.}
  \item{Xr}{A list as \kbd{X}. Each column of a matrix corresponds to the upper interval of the curve effect.}
}

\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}

\seealso{
\code{\link{clustEff}}, for clustering algorithm; \code{\link{summary.clustEff}} and \code{\link{plot.clustEff}}, for summarizing and plotting \code{clustEff} objects.
}
\examples{

# using simulated data

# see the documentation for 'clustEff-package'

}
