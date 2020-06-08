\name{fpcac}
\alias{fpcac}
\title{
Functional Principal Components Analysis Clustering
}
\description{
This function implements the algorithm FPCAC for curves clustering as a variant of a k-means algorithm based on the principal component rotation of data
}
\usage{
fpcac(X, K = 2, fd = NULL, nbasis = 5, norder = 3, nharmonics = 3,
      alpha = 0, niter = 30, Ksteps = 25, conf.level = 0.9, seed, disp = FALSE)
}
\arguments{
\item{X}{Matrix of `curves' of dimension \kbd{n} x \kbd{q}.}
\item{K}{the number of  clusters.}
\item{fd}{If not NULL it overrides X and must be an object of class fd.}
\item{nbasis}{an integer variable specifying the number of basis functions. The default value is 5.}
\item{norder}{an integer specifying the order of b-splines, which is one higher than their degree. The default value is 3.}
\item{nharmonics}{the number of harmonics or principal components to use. The default value is 3.}
\item{alpha}{trimming size, that is  the given proportion of observations to be
discarded.}
\item{niter}{the number or random restarting (larger values provide more accurate solutions.}
\item{Ksteps}{the number of k-mean steps (not too many ksteps are needed).}
\item{conf.level}{the confidence level required.}
\item{seed}{the seed used for reproducibility.}
\item{disp}{if TRUE, it is used to print some information across the algorithm.}
}
\details{FPCAC is a functional
PCA-based clustering approach that provides a variation of the algorithm for curves
clustering proposed by Garcia-Escudero and Gordaliza (2005).

The starting point of the proposed FPCAC is to find a
linear approximation of each curve by a finite $p$ dimensional
vector of coefficients defined by the FPCA scores.

The number of starting clusters \kbd{k} is obtained on the basis of  the scores
volume, such that we assign events to the clusters defined by
events that have a distance less than a fixed threshold (e.g.
90-th percentile) in the space of PCA scores. Once \kbd{k} is obtained
we use a modified version of the trimmed \kbd{k}-means algorithm, that
considers the matrix of FPCA scores instead of the coefficients of
a linear fitting to B-spline bases.

The trimmed \kbd{k}-means clustering algorithm looks for the \kbd{k} centers
\eqn{C_1, ..., C_k} that are solution of the minimization problem:
\deqn{O_k(\alpha)=\min_Y \min_{C_1, \cdots, C_k} \frac{1}{[n(1-\alpha)]}
\sum_{X_i \in Y} \inf_{1\leq j \leq k} || X_i- C_j||^2}

We think that the proposed
approach has the advantage of an immediate use of PCA for
functional data avoiding some objective choices related to spline
fitting as in RCC. Simulations and applications suggest also the
well behavior of the FPCAC algorithm, both in terms of stable and
easily interpretable results.
}
\value{
An object of class \dQuote{\code{fpcac}}, a list containing the following items:
\item{call}{the matched call.}
\item{obj.function}{The percentiles used in the quantile regression coefficient modeling or  objective function \kbd{O_k(\alpha)}.}
\item{centers}{The curves matrix.}
\item{radius}{The vector of clusters.}
\item{clusters}{The mean curves matrix of dimension \kbd{n} x \kbd{k}.}
\item{Xorig}{The atrix of `curves' of dimension \kbd{n} x \kbd{q}.}
\item{fd}{The object obtained by the call of FPCA of class `fd'}
\item{X}{The matrix of `curves' transformed through FPCA of dimension \kbd{p} x \kbd{nharmonics}.}
\item{X.mean}{The mean curves matrix of dimension \kbd{n} x \kbd{k}.}
\item{diss.matrix}{The Euclidean distance matrix of the transformed curves.}
\item{oggSilhouette}{An object of class `silhouette'.}
}
\references{
Adelfio, G., Chiodi, M., D'Alessandro, A. and Luzio, D. (2011) FPCA algorithm for waveform clustering. Journal of Communication and Computer, 8(6), 494-502.

Adelfio, G., Chiodi, M., D'Alessandro, A., Luzio, D., D'Anna, G., Mangano, G. (2012) Simultaneous seismic wave clustering and registration. Computers & Geosciences 44, 60-69.

Garcia-Escudero, L. A. and Gordaliza, A. (2005). A proposal for
robust curve clustering, Journal of classification, 22, 185-201.

}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{opt.fpcac}}.
}
\examples{
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

obj <- fpcac(Y, K = 4, disp = FALSE)
obj
}
\keyword{models}
\keyword{clustering algorithm}
