\name{opt.fpcac}
\alias{opt.fpcac}
\title{
Optimal cluster selection in Functional Principal Components Analysis Clustering
}
\description{
This function provides the optimal selection of clusters for  the algorithm FPCAC,  as a variant of a k-means algorithm based on the principal component rotation of data
}
\usage{
opt.fpcac(X, k.max = 5, method = c("silhouette", "wss"),
          fd = NULL, nbasis = 5, norder = 3, nharmonics = 3,
          alpha = 0, niter = 30, Ksteps = 10, seed,
          diss = NULL, trace=FALSE)
}
\arguments{
\item{X}{Matrix of `curves' of dimension \kbd{n} x \kbd{q}. }
\item{k.max}{the number of cluster used in the optimization step to select the optimal one.}
\item{method}{the method used to select the optimal number of clusters, "silhouette" or "wss" (whithin sum of squares.}
\item{fd}{If not NULL it overrides X and must be an object of class fd.}
\item{nbasis}{an integer variable specifying the number of basis functions. The default value is 5.}
\item{norder}{an integer specifying the order of b-splines, which is one higher than their degree. The default value is 3.}
\item{nharmonics}{the number of harmonics or principal components to use. The default value is 3.}
\item{alpha}{trimming size, that is  the given proportion of observations to be
discarded.}
\item{niter}{the number or random restarting (larger values provide more accurate solutions.}
\item{Ksteps}{the number of k-mean steps (not too many ksteps are needed).}
\item{seed}{the seed used for reproducibility.}
\item{diss}{the dissimilarity matrix used to compute measures "silhouette" or "wss".}
\item{trace}{if TRUE, it is used to print some information across the algorithm.}
}
\details{Silhouette is a method for validate the consistency within clusters, providing a measure  of how similar an object is to its own cluster  compared to other clusters.
The silhouette score S belongs to the interval [-1,1].
 S close to one means that the data is appropriately clustered. If S is close to negative one, datum should be clustered in its neighbouring cluster. S near zero means that the datum is on the border of two natural clusters.

The wss is obtained as the classical sum of the squared deviations from each observation and the cluster centroid, providing  a measure of the variability of the observations within each cluster. Clusters with higher values exhibit greater variability of the observations within the cluster.

}
\value{
a list containing the following items:
\item{obj.function}{the sequence of objective functions.}
\item{clusters}{the matrix in which each columns identify clusters for each fixed K.}
\item{K}{the sequence of K used.}
\item{K.opt}{the optimal number of clusters}
\item{plot}{a ggplot object to plot the curve of silhouette or whithin sum of squares.}
}
\references{ Peter J. Rousseeuw (1987). Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis. Computational and Applied Mathematics. 20, 53-65


 K. V. Mardia, J. T. Kent and J. M. Bibby (1979). Multivariate Analysis. Academic Press.

}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{fpcac}}.
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

num.clust <- opt.fpcac(Y)
obj2 <- fpcac(Y, K = num.clust$K.opt, disp = FALSE)
obj2

}
\keyword{models}
\keyword{clustering algorithm}
