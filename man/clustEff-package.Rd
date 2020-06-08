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
Version: \tab 0.2.0\cr
Date: \tab 2020-05-31\cr
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
\donttest{
clustEff(Beta, Beta.lower = NULL, Beta.upper = NULL,
         k = c(2, min(5, (ncol(Beta)-1))), ask = FALSE, diss.mat, alpha = .5,
         step = c("both", "shape", "distance"),
         cut.method = c("mindist", "length", "conf.int"),
         method = "ward.D2", approx.spline = FALSE, nbasis = 50,
         conf.level = 0.95, stand = FALSE, plot = TRUE, trace = TRUE)

fpcac(X, K = 2, fd = NULL, nbasis = 5, norder = 3, nharmonics = 3,
      alpha = 0, niter = 30, Ksteps = 25, conf.level = 0.95, seed, disp = FALSE)
}

}



