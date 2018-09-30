#' Function for spatial variable selection
#'
#' Perform the variable selection for the spatial Poisson regression model.
#' @param dat.obj List, input data. Must contains:
#' \enumerate{
#'   \item \code{X} numeric matrix, the covariates.
#'   \item \code{y} integer vector, the output count.
#'   \item \code{dist} numeric matrix, the distance matrix.
#'   \item \code{offset} numeric vector, the offset item.
#' }
#' @param alpha.vec numeric vector, vector of possible values of regularization parameter multiplies the adaptive elastic net penalty, valid in [0,1].
#' @param lambda.vec numeric vector, vector of possible values of regularization parameters for the weight L1 norm in the adaptive elastic net penalty, valid in [0,1].
#' @param method string, the method to be used. Options are:
#' \enumerate{
#'   \item \code{"PQL"} penalized quasi-likelihood method consider spatial correlation.
#'   \item \code{"PQL.nocor"} penalized quasi-likelihood method ignore spatial correlation.
#'   \item \code{"APL"} approximate penalized loglikelihood method consider spatial correlation.
#'   \item \code{"APL.nocor"} approximate penalized loglikelihood method ignore spatial correlation.
#' }
#' @param plots bool, if \code{True} plot AIC/BIC parameter contour plots are printed.
#' @param intercept bool, if \code{True} intercept item will be included in model.
#' @param verbose bool, If \code{True}, then various updates are printed during each iteration of the algorithm.
#' @return A list of 13 items:
#' \enumerate{
#'   \item \code{dat.obj} List, a copy of \code{dat.onj} input.
#'   \item \code{start} Initial value given by GLMPQL, transferred to elastic net.
#'   \item \code{L.obj} Linear coefficients' values under each \code{alpha.vec} and \code{lambda.vec}, given by adaptive elastic net.
#'   \item \code{Lout.obj} AIC and BIC values under each \code{L.obj value}, given by adaptive elastic net.
#'   \item \code{contour.out.obj} Data used to plot contour plot with input \code{alpha.vec} and \code{lambda.vec} and output AIC or BIC. Used to choose penalty parameter, given by adaptive elastic net.
#'   \item \code{L.best.obj} Best chosen \code{alpha.vec} and \code{lambda.vec} under \code{contour.out.obj}, given by adaptive elastic net.
#'   \item \code{Lout.best.obj} Best AIC and BIC value under \code{L.best.obj}.
#'   \item \code{L.EN.obj, Lout.EN.obj, contour.out.EN.obj, L.EN.best.obj} similar item but given by elastic net NOT adaptive elastic net, transfer to adaptive elastic net.
#'   \item \code{lasso.weight} Numeric, vector of output adaptive Lasso weight, given by adaptive elastic net.
#'   \item \code{method} String, a copy of the method used for likelihood.
#' }
#' @keywords function
#' @import MASS
#' @import nlme
#' @import fields
#' @import glmnet
#' @import glmmLasso
#' @import stats
#' @import graphics
#' @export
#' @examples
#' #use small.test.dat as the input to fit the spatial Poission regression model
#' test.fit<-SpatialVS(dat.obj=small.test.dat, alpha.vec=0.5,
#' lambda.vec=5, method="PQL", intercept=T, verbose=FALSE)
#####
SpatialVS=function(dat.obj, alpha.vec=seq(0.6, 1, by=0.05), lambda.vec=seq(0.15, 1, len=50), method="PQL", plots=F, intercept=T, verbose=T)
{
  #check for NA's
  
  if(any(is.na(dat.obj$y))|any(is.na(dat.obj$X))|any(is.na(dat.obj$dist))|any(is.na(dat.obj$offset)))
  {
    stop("check for NA's in data object.")
  }

  if(method=="PQL")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testPQL
    SpatialVS.out.FUNCTION=SpatialVS.out
  }

  if(method=="PQL.nocor")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testPQL.nocor
    SpatialVS.out.FUNCTION=SpatialVS.out.nocor
  }

  if(method=="APL")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testLP
    SpatialVS.out.FUNCTION=SpatialVS.out
  }

  if(method=="APL.nocor")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testLP.nocor
    SpatialVS.out.FUNCTION=SpatialVS.out.nocor
  }


  #browser()

  if(verbose)
  {
    cat("Finding initial values \n")
  }


  #browser()

  start.val.ini<-try(spatialVS_simu_ini(simu.data=dat.obj, intercept=intercept, verbose=verbose))


  if(class(start.val.ini)=="try-error")
  {
    res=NULL
    return(res)
  }

  if(verbose)
  {
    cat("ini.beta:", start.val.ini$beta, "ini.theta:", start.val.ini$theta, "\n")
  }

  #lasso.weight=start.val.ini$lasso.weight

  #save(start.val.ini, file="tmp.debug.start.val.ini")
  #load(file="tmp.debug.start.val.ini")

  ###no adpative, i.e., EN

  if(verbose)
  {
    cat("Starting Major looping for alpha and lambda with EN. \n")
  }

  L.EN.obj<-try(SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=0, adaptive=F, intercept=intercept, verbose=verbose), silent=T)

  if(class(L.EN.obj)=="try-error")
  {
    res=NULL
    return(res)
  }

  #browser()
  Lout.EN.obj<-try(SpatialVS.out.FUNCTION(obj=L.EN.obj, data=dat.obj, cc=0), silent=T)

  if(class(Lout.EN.obj)=="try-error")
  {
    res=NULL
    return(res)
  }

  ##########select best tuning parameters
  contour.out.EN.obj=SpatialVS.obj.contour.plot(SpatialVS.obj=L.EN.obj, SpatialVS.out.obj=Lout.EN.obj, type="BIC", plots=plots)

  alpha.EN.best.val=contour.out.EN.obj$alpha.best.val
  lambda.EN.best.val=contour.out.EN.obj$lambda.best.val

  if(verbose)
  {
    cat("Computing best alpha and lambda with EN. \n")
  }

  ###
  L.EN.best.obj=SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.EN.best.val, lambda.vec=lambda.EN.best.val, cc=0, adaptive=F, intercept=intercept, verbose=verbose)

  lasso.weight=as.vector(L.EN.best.obj$beta)
  nn=nrow(dat.obj$X)
  lasso.weight[lasso.weight==0]=(1/nn)

  if(verbose)
  {
    cat("start.ini beta:", round(start.val.ini$beta, 4),"\n")
    cat("EN lasso weight:", round(lasso.weight, 4),"\n")
  }

  if(verbose)
  {
    cat("Starting Major looping for alpha and lambda with AEN. \n")
  }

  L.obj<-try(SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=0, adaptive=T, lasso.weight=lasso.weight, intercept=intercept, verbose=verbose), silent=T)

  if(class(L.obj)=="try-error")
  {
    res=NULL
    return(res)
  }


  if(verbose)
  {
    cat("Computing best alpha and lambda AEN. \n")
  }

  Lout.obj<-try(SpatialVS.out.FUNCTION(obj=L.obj, data=dat.obj, cc=0), silent=T)

  if(class(Lout.obj)=="try-error")
  {
    res=NULL
    return(res)
  }

  ##########select best tuning parameters
  contour.out.obj=SpatialVS.obj.contour.plot(SpatialVS.obj=L.obj, SpatialVS.out.obj=Lout.obj, type="BIC", plots=plots)

  alpha.best.val=contour.out.obj$alpha.best.val
  lambda.best.val=contour.out.obj$lambda.best.val

  #########refit the best model
  L.best.obj=SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.best.val, lambda.vec=lambda.best.val, cc=0, adaptive=T, lasso.weight=lasso.weight, intercept=intercept, verbose=verbose)

  Lout.best.obj=SpatialVS.out(obj=L.best.obj, data=dat.obj, cc=0)

  res=list(dat.obj=dat.obj, start=start.val.ini, L.obj=L.obj, Lout.obj=Lout.obj, contour.out.obj=contour.out.obj, L.best.obj=L.best.obj, Lout.best.obj=Lout.best.obj, L.EN.obj=L.EN.obj, Lout.EN.obj=Lout.EN.obj, contour.out.EN.obj=contour.out.EN.obj, L.EN.best.obj=L.EN.best.obj, lasso.weight=lasso.weight, method=method)

  return(res)
}


