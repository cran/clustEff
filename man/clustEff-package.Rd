\name{clustEff-package}
\alias{clustEff-package}
\docType{package}
\title{
Clusters of effects curves
}
\description{
This package implements a general algorithm to cluster coefficient functions (i.e. clusters of effects) obtained from a quantile regression (qrcm; Frumento and Bottai, 2016). This algorithm is also used for clustering curves observed in time, as in functional data analysis. The objectives of this algorithm vary with the scenario in which it is used, i.e. in the case of a cluster of effects, in a univariate case the objective may be to reduce its dimensionality or in the multivariate case to group similar effects on a covariate. In the case of a functional data analysis the main objective is to cluster waves or any other function of time or space. Sottile G. and Adelfio G. (2019) <https://doi.org/10.1007/s00180-018-0817-8>.
}
\details{
\tabular{ll}{
Package: \tab clustEff\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2024-01-22\cr
License: \tab GPL-2\cr
}
The function \code{\link{clustEff}} allows to specify the type of the curves to apply the proposed clustering algorithm. The function \code{\link{extract.object}}  extracts the matrices, in case of multivariate response, through the quantile regression coefficient modeling, useful to run the main algorithm.  The auxiliary functions \code{\link{summary.clustEff}} and \code{\link{plot.clustEff}} can be used to extract information from the main algorithm. In the new version of the package you can also find a PCA-based clustering approach called Functional Principal Components Analysis Clustering (FPCAC). Main function of this algorithm is \code{\link{fpcac}}, and some auxiliary functions are \code{\link{summary.fpcac}} and \code{\link{plot.fpcac}}.
}
\author{
Gianluca Sottile

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
}
\references{
Sottile, G., Adelfio, G. \emph{Clusters of effects curves in quantile regression models}. Comput Stat 34, 551–569 (2019). https://doi.org/10.1007/s00180-018-0817-8

Sottile, G and Adelfio, G (2017). \emph{Clustering of effects through quantile regression}. Proceedings 32nd International Workshop of Statistical Modeling, Groningen (NL), vol.2 127-130, https://iwsm2017.webhosting.rug.nl/IWSM_2017_V2.pdf.

Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, doi: 10.1111/biom.12410.

Adelfio, G., Chiodi, M., D'Alessandro, A. and Luzio, D. (2011) \emph{FPCA algorithm for waveform clustering}. Journal of Communication and Computer, 8(6), 494-502.

}
\keyword{ package }
\examples{

# Main functions:
set.seed(1234)
n <- 300
x <- 1:n/n

Y <- matrix(0, n, 30)

sigma2 <- 4*pmax(x-.2, 0) - 8*pmax(x-.5, 0) + 4*pmax(x-.8, 0)

mu <- sin(3*pi*x)
for(i in 1:10) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

mu <- cos(3*pi*x)
for(i in 11:23) Y[,i] <- mu + rnorm(length(x), 0, pmax(sigma2,0))

mu <- sin(3*pi*x)*cos(pi*x)
for(i in 24:28) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

mu <- 0 #sin(1/3*pi*x)*cos(2*pi*x)
for(i in 29:30) Y[, i] <- mu + rnorm(length(x), 0, pmax(sigma2, 0))

clustEff(Y)

fpcac(Y, K = opt.fpcac(Y)$K.opt)
}



