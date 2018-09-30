\name{SpatialVS}
\alias{SpatialVS}
\title{Function for spatial variable selection}
\usage{
SpatialVS(dat.obj, alpha.vec = seq(0.6, 1, by = 0.05),
  lambda.vec = seq(0.15, 1, len = 50), method = "PQL", plots = F,
  intercept = T, verbose = T)
}
\arguments{
\item{dat.obj}{List, input data. Must contains:
\enumerate{
  \item \code{X} numeric matrix, the covariates.
  \item \code{y} integer vector, the response in counts.
  \item \code{dist} numeric matrix, the distance matrix.
  \item \code{offset} numeric vector, the offset item.
}}

\item{alpha.vec}{numeric vector, a vector of possible values of regularization parameter. The range is [0,1].}

\item{lambda.vec}{numeric vector, a vector of positive values of regularization parameter. }

\item{method}{string, the method to be used. Options are:
\enumerate{
  \item \code{"PQL"} penalized quasi-likelihood method that considers spatial correlation.
  \item \code{"PQL.nocor"} penalized quasi-likelihood method that ignores spatial correlation.
  \item \code{"APL"} approximate penalized loglikelihood method that considers spatial correlation.
  \item \code{"APL.nocor"} approximate penalized loglikelihood method that ignores spatial correlation.
}}

\item{plots}{bool, if \code{True}, contour plot of AIC/BIC values is generated.}

\item{intercept}{bool, if \code{True}, an intercept item will be included in model.}

\item{verbose}{bool, if \code{True}, various updates are printed during each iteration of the algorithm.}
}
\value{
A list of 13 items:
\enumerate{
  \item \code{dat.obj}, List, a copy of the \code{dat.obj} input.
  \item \code{start}, Initial values of parameters given by glmmPQL().
  \item \code{L.obj},  Regression coefficients under each \code{alpha.vec} and \code{lambda.vec}, under the adaptive elastic net.
  \item \code{Lout.obj}, AIC and BIC values under each \code{L.obj value}, under the adaptive elastic net.
  \item \code{contour.out.obj}, Object used to generate contour plot as a function of \code{alpha.vec} and \code{lambda.vec}, with  AIC or BIC as the output. Used to choose best penalty parameter, under the adaptive elastic net.
  \item \code{L.best.obj}, Model fitting results under the best chosen \code{alpha.vec} and \code{lambda.vec}, under the adaptive elastic net.
  \item \code{Lout.best.obj}, Best BIC value for \code{L.best.obj}.
  \item \code{L.EN.obj, Lout.EN.obj, contour.out.EN.obj, L.EN.best.obj}, Similar items but under the elastic penalty.
  \item \code{lasso.weight}, Numeric, specifies the adaptive Lasso weight.
  \item \code{method}, String, the method used for computing the approximate likelihood function.
}
}
\description{
Perform variable selection for the spatial Poisson regression model under adaptive elastic net penalty.
}

\references{
Xie, Y., Xu, L., Li, J., Deng, X., Hong, Y., Kolivras, K., and Gaines, D. N. (2018). Spatial Variable Selection and An Application to Virginia Lyme Disease Emergence. Preprint, arXiv:1809.06418 [stat.AP].
}

\examples{
#use small.test.dat as the input to fit the spatial Poission regression model.
#a grid of alpha.vec and lambda.vec is typically used.
#Here one point of alpha.vec and lambda.vec is used for fast illustration.

test.fit<-SpatialVS(dat.obj=small.test.dat, alpha.vec=0.5,
lambda.vec=5, method="PQL", intercept=TRUE, verbose=FALSE)
}
\keyword{function}
