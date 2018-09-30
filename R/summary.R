#' Function for spatial variable selection's summary
#'
#' return the summary of a \link{SpatialVS}'s return
#'
#' @param obj List, returned by \link{SpatialVS}.
#'
#' @return A detailed summary of \link{SpatialVS} object.
#' @export
#' @examples
#' library(MASS)
#' library(nlme)
#' test.fit<-SpatialVS(dat.obj=small.test.dat, alpha.vec=0.5, lambda.vec=5,
#' method="PQL", intercept=T, verbose=FALSE)
#' SpatialVS.summary(test.fit)


SpatialVS.summary=function(obj)
{
  dat.obj=obj$data.obj
  L.best.obj=obj$L.best.obj
  Lout.best.obj=obj$Lout.best.obj

  res=c(L.best.obj$beta.res, L.best.obj$theta.res)

  return(res)
}
