\name{clustEff}
\alias{clustEff}
\title{
Cluster Effects Algorithm
}
\description{
This function implements the algorithm to cluster curves of effects obtained from a
quantile regression (qrcm; Frumento and Bottai, 2015) in which the coefficients are described by flexible parametric functions of the order of the quantile. This algorithm can be also used for clustering of curves observed in time, as in functional data analysis.
}
\usage{
clustEff(Beta, Beta.lower = NULL, Beta.upper = NULL,
         k = c(2, min(5, (ncol(Beta)-1))), ask = FALSE, diss.mat, alpha = .5,
         step = c("both", "shape", "distance"),
         cut.method = c("mindist", "length", "conf.int"),
         method = "ward.D2", approx.spline = FALSE, nbasis = 50,
         conf.level = 0.9, stand = FALSE, plot = TRUE, trace = TRUE)
}
\arguments{
\item{Beta}{A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of curves to cluster and \kbd{n} is either the length of percentiles used in the quantile regression or the length of the time vector.}
\item{Beta.lower}{A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of lower interval of the   curves to cluster and \kbd{n} the length of percentiles used in quantile regression. Used only if cluster.effects=TRUE.}
\item{Beta.upper}{A matrix \kbd{n} x \kbd{q}. \kbd{q} represents the number of upper interval of the curves to cluster and \kbd{n} the length of percentiles used in quantile regression. Used only if cluster.effects=TRUE.}
\item{k}{It represents the number of clusters to look for. If it is two-length vector (k.min - k.max) an optimization is performed, if it is a unique value it is fixed.}
\item{ask}{If TRUE, after plotting the dendrogram, the user make is own choice about how many cluster to use.}
\item{diss.mat}{a dissimilarity matrix, obtained by using distshape function.}
\item{alpha}{It is the alpha-percentile used for computing the dissimilarity matrix. The default value is alpha=.5.}
\item{step}{The steps used in computing the dissimilarity matrix. Default is "both"=("shape" and  "distance")}
\item{cut.method}{The method used in optimization step to look for the optimal number of clusters. Default is "mindist", however if Beta.lower and Beta.upper are available the suggested method is "conf.int".}
\item{method}{The agglomeration method to be used.}
\item{approx.spline}{If TRUE, Beta is approximated by a smooth spline.}
\item{nbasis}{An integer variable specifying the number of basis functions. Only when approx.spline=TRUE}
\item{conf.level}{the confidence level required.}
\item{stand}{If TRUE, the argument Beta is standardized.}
\item{plot}{If TRUE, dendrogram, boxplot and clusters are plotted.}
\item{trace}{If TRUE, some informations are printed.}

}
\details{
  Quantile regression models conditional quantiles of a response variabile,
  given a set of covariates. Assume that each coefficient can be expressed as a parametric function of      \eqn{p} in the form: \deqn{\beta(p | \theta) = \theta_{0} + \theta_1 b_1(p) + \theta_2 b_2(p) +      \ldots}{ \beta(p | \theta) = \theta0 + \theta1*b1(p) + \theta2*b2(p) + \ldots}
  where \eqn{b_1(p), b_2(p, \ldots)}{b1(p), b2(p), \ldots} are known functions of \eqn{p}.
}
\value{
An object of class \dQuote{\code{clustEff}}, a list containing the following items:
\item{call}{the matched call.}
\item{p}{The percentiles used in quantile regression coefficient modeling or the time otherwise.}
\item{X}{The curves matrix.}
\item{clusters}{The vector of clusters.}
\item{X.mean}{The mean curves matrix of dimension \kbd{n} x \kbd{k}.}
\item{X.mean.dist}{The within cluster distance from the mean curve.}
\item{X.lower}{The lower bound matrix.}
\item{X.mean.lower}{The mean lower bound of dimension \kbd{n} x \kbd{k}.}
\item{X.upper}{The upper bound matrix.}
\item{X.mean.upper}{The mean upper bound of dimension \kbd{n} x \kbd{k}.}
\item{Signif.interval}{The matrix of dimension \kbd{n} x \kbd{k} containing the intervals in which each mean lower and upper bounds don't include the zero.}
\item{k}{The number of selected clusters.}
\item{diss.matrix}{The dissimilarity matrix.}
\item{X.mean.diss}{The within cluster dissimilarity.}
\item{oggSilhouette}{An object of class \dQuote{\code{silhouette}}.}
\item{oggHclust}{An object of class \dQuote{\code{hclust}}.}
\item{distance}{A vector of goodness measures used to select the best number of clusters.}
\item{step}{The selected step.}
\item{method}{The used agglomeration method.}
\item{cut.method}{The used method to select the best number of clusters.}
\item{alpha}{The selected alpha-percentile.}
}
\references{
Sottile, G., Adelfio, G. Clusters of effects curves in quantile regression models. Comput Stat 34, 551–569 (2019). https://doi.org/10.1007/s00180-018-0817-8

Sottile, G and Adelfio, G (2017). \emph{Clustering of effects through quantile regression}.
Proceedings 32nd International Workshop of Statistical Modeling, Groningen (NL), vol.2 127-130, https://iwsm2017.webhosting.rug.nl/IWSM_2017_V2.pdf.

Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, doi: 10.1111/biom.12410.
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
\seealso{
\code{\link{summary.clustEff}}, \code{\link{plot.clustEff}},
for summary and plotting.
\code{\link{extract.object}} to extract useful objects for the clustering algorithm through a quantile regression coefficient modeling in a multivariate case.
}
\examples{

# CURVES EFFECTS CLUSTERING

set.seed(1234)
n <- 300
q <- 2
k <- 5
x1 <- runif(n, 0, 5)
x2 <- runif(n, 0, 5)

X <- cbind(x1, x2)
rownames(X) <- 1:n
colnames(X) <- paste0("X", 1:q)

theta1 <- matrix(c(1, 1, 0, 0, 0, .5, 0, .5, 1, 2, .5, 0, 2, 1, .5),
                 ncol=k, byrow=TRUE)

theta2 <- matrix(c(1, 1, 0, 0, 0, -.3, 0, .5, 1, .5, -1.5, 0, -1, -.5, 1),
                 ncol=k, byrow=TRUE)

theta3 <- matrix(c(1, 1, 0, 0, 0, .3, 0, -.5, -1, 2, -.5, 0, 1, -.5, -1),
                 ncol=k, byrow=TRUE)

rownames(theta3) <- rownames(theta2) <- rownames(theta1) <-
    c("(intercept)", paste("X", 1:q, sep=""))
colnames(theta3) <- colnames(theta2) <- colnames(theta1) <-
    c("(intercept)", "qnorm(p)", "p", "p^2", "p^3")

Theta <- list(theta1, theta2, theta3)

B <- function(p, k){matrix(cbind(1, qnorm(p), p, p^2, p^3), nrow=k, byrow=TRUE)}
Q <- function(p, theta, B, k, X){rowSums(X * t(theta \%*\% B(p, k)))}

Y <- matrix(NA, nrow(X), 15)
for(i in 1:15){
  if(i <= 5) Y[, i] <- Q(runif(n), Theta[[1]], B, k, cbind(1, X))
  if(i <= 10 & i > 5) Y[, i] <- Q(runif(n), Theta[[2]], B, k, cbind(1, X))
  if(i <= 15 & i > 10) Y[, i] <- Q(runif(n), Theta[[3]], B, k, cbind(1, X))
}

XX <- extract.object(Y, X, intercept=TRUE, formula.p= ~ I(p) + I(p^2) + I(p^3))

obj <- clustEff(XX$X$X1, Beta.lower=XX$Xl$X1, Beta.upper=XX$Xr$X1, cut.method = "conf.int")
summary(obj)
plot(obj, xvar="clusters", col = 1:3)
plot(obj, xvar="dendrogram")
plot(obj, xvar="boxplot")

obj2 <- clustEff(XX$X$X2, Beta.lower=XX$Xl$X2, Beta.upper=XX$Xr$X2, cut.method = "conf.int")
summary(obj2)
plot(obj2, xvar="clusters", col=1:3)
plot(obj2, xvar="dendrogram")
plot(obj2, xvar="boxplot")


\dontrun{
set.seed(1234)
n <- 300
q <- 15
k <- 5
X <- matrix(rnorm(n*q), n, q); X <- scale(X)
rownames(X) <- 1:n
colnames(X) <- paste0("X", 1:q)

Theta <- matrix(c(1, 1, 0, 0, 0,
                  .5, 0, .5, 1, 1,
                  .5, 0, 1, 2, .5,
                   .5, 0, 1, 1, .5,
                  .5, 0, .5, 1, 1,
                   .5, 0, .5, 1, .5,
                 -1.5, 0, -.5, 1, 1,
                  -1, 0, .5, -1, -1,
                 -.5, 0, -.5, -1, .5,
                  -1, 0, .5, -1, -.5,
                -1.5, 0, -.5, -1, -.5,
                  2, 0, 1, 1.5, 2,
                  2, 0, .5, 1.5, 2,
                  2.5, 0, 1, 1, 2,
                  1.5, 0, 1.5, 1, 2,
                  3, 0, 2, 1, .5),
                 ncol=k, byrow=TRUE)
rownames(Theta) <- c("(intercept)", paste("X", 1:q, sep=""))
colnames(Theta) <- c("(intercept)", "qnorm(p)", "p", "p^2", "p^3")

B <- function(p, k){matrix(cbind(1, qnorm(p), p, p^2, p^3), nrow=k, byrow=TRUE)}
Q <- function(p, theta, B, k, X){rowSums(X * t(theta \%*\% B(p, k)))}

s <- matrix(1, q+1, k)
s[2:(q+1), 2] <- 0
s[1, 3:k] <- 0

Y <- Q(runif(n), Theta, B, k, cbind(1, X))
XX <- extract.object(Y, X, intercept = TRUE, formula.p= ~ I(p) + I(p^2) + I(p^3))

obj3 <- clustEff(XX$X, Beta.lower=XX$Xl, Beta.upper=XX$Xr, cut.method = "conf.int")
summary(obj3)

# changing the alpha-percentile clusters are correctly identified

obj4 <- clustEff(XX$X, Beta.lower=XX$Xl, Beta.upper=XX$Xr, cut.method = "conf.int",
                 alpha = 0.25)
summary(obj4)

# CURVES CLUSTERING IN FUNCTIONAL DATA ANALYSIS

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

obj5 <- clustEff(Y)
summary(obj5)
plot(obj5, xvar="clusters", col=1:4)
plot(obj5, xvar="dendrogram")
plot(obj5, xvar="boxplot")
}

}
\keyword{models}
\keyword{clustering algorithm}
