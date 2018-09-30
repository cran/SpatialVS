#Performs penalized quasi-likelihood method that considers spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS

SpatialVS_P_testPQL=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist

  n=length(y)

  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_PQLpenb
         })

  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {

      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }

      nIter=0

      beta=start$beta
      theta=start$theta
      D=covStruct(theta, cc=cc, dist=dist)
      convCov=convPar=1



      p=length(beta)

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }


      b=start$b

      ###main event
      while ((control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {


        #cat("nIter=", nIter, "convCov=", convCov, "\n")

        nIter_beta=1
        convPar=1

        while((nIter_beta<2) & (convPar > control$tol2))   #nIter_beta<30
        {
          ###first step: get the current estimates of beta and b using iwls
          res1=iwls(X,y,D,offset,beta,b,tol=control$iwls,dist=dist)

          ####second step: the elastic-net step
          beta.new=res1$beta;
          mu.new=res1$mu;
          b.new=res1$b

          #print(mean(b.new))

          Z=X%*%beta.new-rep(1,n)+y/mu.new

          if(adaptive)
          {

            pf1.wts=1/(abs(lasso.weight))
            pf2.wts=rep(1,p)

            if(intercept)
            {
              pf1.wts[1]=0
              pf2.wts[1]=0
            }

            res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)
            beta.new=res2$coef

          }else{
                  pf1.wts = rep(1,p)
                  pf2.wts=rep(1,p)

                  if(intercept)
                  {
                     pf1.wts[1]=0
                     pf2.wts[1]=0
                  }

                  res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)

                  beta.new=res2$coef

                }


          convPar <- sum(sqrt(crossprod(beta.new-beta))/(1+sqrt(crossprod(beta.new))))

          #cat("        nIter_beta=", nIter_beta, "convPar=", convPar, "\n")

          nIter_beta=nIter_beta+1

          #print(convPar)
          #print(beta.new)

          beta=beta.new
          b=b.new


        }#end of beta iteration

        ###third step: update the theta

        mu=as.vector(exp(X%*%beta + b + offset))
        ystar=as.vector(X%*%beta+b+(y-mu)/mu)


        nonzerobeta=beta[beta!=0]
        X_L=as.matrix(X[, beta!=0])

        if(adaptive)
        {

          nonzerolasso.weight=abs(lasso.weight)[beta!=0]

          Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta)*abs(nonzerolasso.weight))+2*lambda.vec[j]*(1-alpha.vec[i])


          if(intercept & beta[1]!=0)
          {
             Sigma_diag[1]=0
          }



        }else{
               Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta))+2*lambda.vec[j]*(1-alpha.vec[i])

               if(intercept & beta[1]!=0)
               {
                 Sigma_diag[1]=0
               }


             }

        Sigma_pen=diag(Sigma_diag, nrow=length(Sigma_diag), ncol=length(Sigma_diag))


        system.time(res3<-optim(theta, cov.fun, mu=mu, nonzerobeta=nonzerobeta, n=n, X_L=X_L, Sigma_pen=Sigma_pen,  dist=dist, lower=-5, upper=10, method="L-BFGS-B", cc=cc, offset=offset, X=X, y=y, b.start=b, tol=control.default$iwls, beta=beta))

        theta.new=res3$par

        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        #print(convCov)
        #print(theta.new)

        nIter=nIter+1

        #print(paste("Number of iterations:", nIter))

        ###update
        theta=theta.new;
        D=covStruct(theta,cc=cc,dist=dist)

      }#end of main iteration


      res_b=iwls_b(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
      b=res_b$b

      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])

      if(verbose)
      {
       print(beta)
       print(theta)
       print(nIter)
      }

    }
  }

  ori.theta.res=cbind(exp(theta.res[,1]), exp(theta.res[,2]))

  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))

}

#Performs penalized quasi-likelihood method that ignores spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS

SpatialVS_P_testLP.nocor=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist



  n=length(y)

  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_LP.nocor
         })

  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {
      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }

      nIter=0

      beta=start$beta
      theta=start$theta[1]
      D=covStruct.nocor(theta,n=n)
      convCov=convPar=1
      p=length(beta)

      b=start$b

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }

      pf1.wts=1/abs(lasso.weight)

      if(intercept)
      {
        pf1.wts[1]=0
      }

      ###main event
      while((control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {

        #if(nIter==5) {browser()}


        # iterate beta given theta
        # estimate beta given theta
        nIter_beta=1
        convPar=1

        # estimate beta given theta
        while( (convPar > control$tol2)&(20 > nIter_beta) )
        {
          nIter_beta=nIter_beta+1
          beta.ini=beta

          for(s in 1:p)
          {
            e_s=rep(0, p)
            e_s[s]=1

            res1=iwls_b_new(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)

            beta=res1$beta
            mu=res1$mu
            b=res1$b
            eta=log(mu)

            x_s=as.numeric(X[,s])
            DW=D%*%diag(mu)

            if(intercept & s==1)
            {
               beta_sm=sum(x_s*(mu-y))+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))
               h_sm=sum(x_s*x_s*mu)

            }else{
                   beta_sm=sum(x_s*(mu-y))+2*lambda.vec[j]*(1-alpha.vec[i])*beta[s]+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))
                   h_sm=sum(x_s*x_s*mu)+2*lambda.vec[j]*(1-alpha.vec[i])
                 }

            if(intercept & s==1)
            {
               d_sm=-beta_sm/h_sm

            }else{
                    w_s=pf1.wts[s]
                    d_sm=median(c((w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm, -beta[s], (-w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm))
                 }


            beta=beta+d_sm*e_s

          }
          beta.new=beta
          convPar <- sum(sqrt(crossprod(beta.new-beta.ini))/(1+sqrt(crossprod(beta.new))))
          #print(convPar)

        }



        # estimate covariance parameters
        res3=optim(theta, cov.fun, beta=as.vector(beta.new), X=X, y=y, cc=cc, b.start=as.vector(b), offset=offset,  dist=dist, lower=-5, upper=10, method="L-BFGS-B")

        theta.new=res3$par

        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        #print(convCov)


        nIter=nIter+1

        #print(paste("Number of iterations:", nIter))

        ###update
        beta=beta.new
        theta=theta.new
        D=covStruct.nocor(theta,n=n)

      }

      #print(paste("Number of iterations:", nIter))

      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])

      if(verbose)
      {
        print(beta)
        print(theta)
        print(nIter)
      }

    }
  }

  ori.theta.res=cbind(exp(theta.res))

  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))

}

#Performs pproximate penalized loglikelihood method that considers spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS

SpatialVS_P_testLP=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist



  n=length(y)

  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_LP
         })

  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {

      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }

      nIter=0

      beta=start$beta
      theta=start$theta
      D=covStruct(theta,cc=cc, dist=dist)
      convCov=convPar=1
      p=length(beta)

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }

      pf1.wts=1/abs(lasso.weight)

      if(intercept)
      {
        pf1.wts[1]=0
      }

      b=start$b

      ###main event
      while ( (control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2))
      {



        nIter_beta=1
        convPar=1

        # estimate beta given theta
        while( (convPar > control$tol2)&(20 > nIter_beta) )
        {
          nIter_beta=nIter_beta+1

          beta.ini=beta

          for(s in 1:p)
          {
            e_s=rep(0, p)
            e_s[s]=1

            res1=iwls_b_new(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)

            beta=res1$beta
            mu=res1$mu
            b=res1$b
            eta=log(mu)

            x_s=as.numeric(X[,s])
            DW=D%*%diag(mu)

            if(intercept & s==1)
            {
               beta_sm=sum(x_s*(mu-y))+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))
               h_sm=sum(x_s*x_s*mu)

            }else{
                   beta_sm=sum(x_s*(mu-y))+2*lambda.vec[j]*(1-alpha.vec[i])*beta[s]+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))
                   h_sm=sum(x_s*x_s*mu)+2*lambda.vec[j]*(1-alpha.vec[i])
                 }

            if(intercept & s==1)
            {
               d_sm=-beta_sm/h_sm


            }else{
                    w_s=pf1.wts[s]
                    d_sm=median(c((w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm, -beta[s], (-w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm))
                 }

            beta=beta+d_sm*e_s



          }

          beta.new=beta
          convPar <- sum(sqrt(crossprod(beta.new-beta.ini))/(1+sqrt(crossprod(beta.new))))
          #print(convPar)

        }

        # estimate covariance parameters
        res3=optim(theta, cov.fun, beta=beta.new, X=X, y=y, cc=cc, b.start=b, offset=offset,  dist=dist, lower=-5, upper=10, method="L-BFGS-B")

        theta.new=res3$par

        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        #print(convCov)

        nIter=nIter+1

        #print(paste("Number of iterations:", nIter))

        ###update
        beta=beta.new
        theta=theta.new
        D=covStruct(theta,cc=cc,dist=dist)

      }

      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])

      if(verbose)
      {
        print(beta)
        print(theta)
        print(nIter)
      }

    }

  }

  ori.theta.res=cbind(exp(theta.res[,1]), exp(theta.res[,2]))

  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))

}


#Performs pproximate penalized loglikelihood method that ignores spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS

SpatialVS_P_testPQL.nocor=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist

  n=length(y)

  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_PQLpen.nocor
         })

  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {
      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }


      nIter=0

      beta=start$beta
      theta=start$theta[1]
      D=covStruct.nocor(theta, n=n)
      convCov=convPar=1
      p=length(beta)

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }

      b=start$b

      ###main event
      while ( (control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {


        convPar=1
        while(convPar > control$tol2)
        {
          ###first step: get the current estimates of beta and b using iwls
          res1=iwls(X,y,D,offset,beta,b,tol=control$iwls,dist=dist)

          ####second step: the elastic-net step
          beta.new=res1$beta
          mu.new=res1$mu
          b.new=res1$b

          Z=X%*%beta.new-rep(1,n)+y/mu.new

          if(adaptive)
          {
            # adaptive elastic-net

            pf1.wts=1/(abs(lasso.weight))
            pf2.wts=rep(1,p)

            if(intercept)
            {
              pf1.wts[1]=0
              pf2.wts[1]=0
            }

            res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)

            beta.new=res2$coef

          }else{

                 pf1.wts = rep(1,p)
                 pf2.wts=rep(1,p)

                 if(intercept)
                 {
                   pf1.wts[1]=0
                   pf2.wts[1]=0
                 }

                 res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)

                  beta.new=res2$coef

          }


          convPar <- sum(sqrt(crossprod(beta.new-beta))/(1+sqrt(crossprod(beta.new))))

          beta=beta.new
          b=b.new
          ###third step: update the theta
        }


        mu=as.vector(exp(X%*%beta + b + offset))
        ystar=as.vector(X%*%beta+b+(y-mu)/mu)

        nonzerobeta=beta[beta!=0]
        X_L=X[, beta!=0]

        if(adaptive)
        {
          nonzerolasso.weight=abs(lasso.weight)[beta!=0]
          Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta)*abs(nonzerolasso.weight))+2*lambda.vec[j]*(1-alpha.vec[i])

          if(intercept & beta[1]!=0)
          {
             Sigma_diag[1]=0
          }

        }else{
               Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta))+2*lambda.vec[j]*(1-alpha.vec[i])

               if(intercept & beta[1]!=0)
               {
                 Sigma_diag[1]=0
               }

             }

        Sigma_pen=diag(Sigma_diag, nrow=length(Sigma_diag), ncol=length(Sigma_diag))

        res3=optim(theta, cov.fun, mu=mu, nonzerobeta=nonzerobeta, n=n, X_L=X_L, Sigma_pen=Sigma_pen, ystar=ystar, dist=dist, lower=-5, upper=10, method="L-BFGS-B", cc=cc)

        theta.new=res3$par

        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        nIter=nIter+1

        #print(paste("Number of iterations:", nIter))

        ###update
        theta=theta.new;
        D=covStruct.nocor(theta,n=n)
      }


      res_b=iwls_b(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
      b=res_b$b

      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])

      if(verbose)
      {
        print(beta)
        print(theta)
      }

    }
  }

  ori.theta.res=exp(theta.res)
  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res,
              alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))
}


#Function used to compute the AIC/BIC of spatial fitting from functions such SpatialVS_P_testPQL
#obj is the output of SpatialVS_P_testPQL
#data is the data matrix
#cc is a variable used in computing covariance matrix

SpatialVS.out=function(obj, data, cc)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist

  n=length(y)

  lld=bic=aic=dev=ERIC=rep(0, dim(obj$beta.res)[1])

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  for(i in 1:length(dev))
  {

    tt=as.numeric(X%*%obj$beta.res[i,]) + as.numeric(obj$b.res[i,]) + offset
    D=covStruct(obj$theta.res[i,], cc=cc, dist=dist)

    detextraterm=determinant(D%*%diag(as.vector(exp(tt)))+diag(rep(1,n)))

    if(sum(log(factorial(y)))==Inf)
    {
      lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt) ) -
        as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))
        -0.5*(detextraterm$modulus)*(detextraterm$sign)
    }else{
      lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt) ) - sum(log(factorial(y))) -
        as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))
        -0.5*(detextraterm$modulus)*(detextraterm$sign)
    }

    k=sum(obj$beta.res[i,]!=0)+2

    bic[i]=-2*lld[i]+k*log(n)

    aic[i]=-2*lld[i]+k*2

    dev[i]=sum(2*(y[y!=0]*log(y[y!=0]/exp(tt[y!=0]))))+sum(2*(-(y-exp(tt))))



    ERIC[i]=-2*lld[i]+2*k*log(n/((obj$alpha.res[,i])*(obj$lambda.res[,i])))
  }

  return(list(lld=lld, bic=bic, dev=dev, aic=aic, ERIC=ERIC))

}


#Function used to compute the AIC/BIC of spatial fitting from functions such SpatialVS_P_testPQL.nocor
#obj is the output of SpatialVS_P_testPQL
#data is the data matrix
#cc is a variable used in computing covariance matrix

SpatialVS.out.nocor=function(obj, data, cc)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist

  n=length(y)

  lld=bic=aic=dev=rep(0, dim(obj$beta.res)[1])

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  for(i in 1:length(dev))
  {
    tt=as.numeric(X%*%obj$beta.res[i,]) + as.numeric(obj$b.res[i,]) + offset
    D=covStruct.nocor(obj$theta.res[i,],n=n)

    detextraterm=determinant(D%*%diag(as.vector(exp(tt)))+diag(rep(1,n)))

    lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt-offset) ) -
      as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))-
      0.5*(detextraterm$modulus)*(detextraterm$sign)

    k=sum(obj$beta.res[i,]!=0)+1

    bic[i]=-2*lld[i]+k*log(n)

    aic[i]=-2*lld[i]+k*2

    dev[i]=sum(2*(y[y!=0]*log(y[y!=0]/exp(tt[y!=0]))))+sum(2*(-(y-exp(tt))))
  }

  return(list(lld=lld, bic=bic, dev=dev, aic=aic))

}


#make a contour plot for the AIC or BIC for penalty tuning parameters.

SpatialVS.obj.contour.plot=function(SpatialVS.obj, SpatialVS.out.obj, type="BIC", nlevels=5, sep=5, plots=T)
{

  x=SpatialVS.obj$alpha.vec
  y=SpatialVS.obj$lambda.vec



  switch(type,
         "BIC"={
           z=matrix(SpatialVS.out.obj$bic, nrow=length(x), ncol=length(y), byrow=T)
         },
         "AIC"={
           z=matrix(SpatialVS.out.obj$aic, nrow=length(x), ncol=length(y), byrow=T)
         },
         "Deviance"={
           z=matrix(SpatialVS.out.obj$dev, nrow=length(x), ncol=length(y), byrow=T)
         },
         "LogLikelihood"={
           z=matrix(SpatialVS.out.obj$lld, nrow=length(x), ncol=length(y), byrow=T)
         })

  id.x=order(x)
  id.y=order(y)


  tmp.idx.mat=(z==min(z))
  alpha.best.val=x[rowSums(tmp.idx.mat)>0]
  lambda.best.val=y[colSums(tmp.idx.mat)>0]
  alpha.best.val=alpha.best.val[1]
  lambda.best.val=lambda.best.val[1]

  if(plots)
  {
    contour(x=x[id.x], y=y[id.y], z=z[id.x, id.y],  xlab=expression(alpha), ylab=expression(lambda),levels=seq(min(z)+0.01, min(z)+3,,5), main=type)

    points(alpha.best.val, lambda.best.val, pch=16, col=2, lwd=2)


  }

  #print(c(alpha.best.val, lambda.best.val))

  res=list(x=x, y=y, z=z, alpha.best.val=alpha.best.val, lambda.best.val=lambda.best.val)

  return(invisible(res))
}

#Used by SpatialVS in the estimation of covariance paramters under penalty

covobj_PQLpenb=function(theta, cc, mu, nonzerobeta, beta, b.start, tol=control.default$iwls, n, X_L, X, y, Sigma_pen,  dist, offset)
{
  D=covStruct(theta,cc=cc,dist=dist)

  res1=iwls_b(X=X, y=y, D=D, offset=offset, beta=beta, b=b.start, tol=tol,dist=dist)
  b=res1$b
  mu=res1$mu


  V=diag(1/mu)+D
  solveV=solve(V)
  detV=determinant(V)
  detXVX=determinant(t(X_L)%*%solveV%*%X_L+Sigma_pen)

  ystar=X%*%beta+b+(y-mu)/mu

  res=as.numeric(0.5*detV$modulus*detV$sign+0.5*detXVX$modulus*detXVX$sign+0.5*t(ystar-X_L%*%nonzerobeta)%*%solveV%*%(ystar-X_L%*%nonzerobeta))

  return(res)
}

#Used by SpatialVS in the estimation of covariance paramters under penalty when ignores spatial correlation

covobj_PQLpen.nocor=function(theta, cc, mu, nonzerobeta, n, X_L, Sigma_pen, ystar, dist)
{

  V=diag((1/mu)+exp(theta))

  solveV=diag(1/((1/mu)+exp(theta)))

  detV=sum(log((1/mu)+exp(theta)))
  detXVX=determinant(t(X_L)%*%solveV%*%X_L+Sigma_pen)

  return(as.numeric(0.5*detV+0.5*detXVX$modulus*detXVX$sign+
                      0.5*t(ystar-X_L%*%nonzerobeta)%*%solveV%*%(ystar-X_L%*%nonzerobeta)))
}


#Used by SpatialVS in the estimation of covariance paramters under spatial correlation

covobj_LP=function(theta, cc, dist, X, y, offset=offset, beta, b.start, tol=control.default$iwls)
{
  D=covStruct(theta,cc=cc,dist=dist)
  n=length(y)

  b=b.start
  mu=exp(X%*%as.vector(beta)+b+offset)
  mu=as.vector(mu)


  detextraterm=determinant(D%*%diag(as.vector(mu))+diag(rep(1,n)))

  res=sum(mu)-sum(y*(X%*%as.vector(beta)+b))+0.5*t(b)%*%solve(D)%*%b+0.5*(detextraterm$modulus)*(detextraterm$sign)

  return(as.numeric(res))

}


#Used by SpatialVS in the estimation of covariance paramters under no spatial correlation

covobj_LP.nocor=function(theta, cc, dist, X, y, offset=offset, beta, b.start, tol=control.default$iwls)
{
  n=length(y)

  D=covStruct.nocor(theta,n=n)

  b=b.start
  eta=X%*%as.vector(beta)+b+offset
  eta=as.vector(eta)
  mu=exp(eta)

  detextraterm=determinant(D%*%diag(as.vector(mu))+diag(rep(1,n)))

  res=0.5*t(b)%*%solve(D)%*%b +0.5*(detextraterm$modulus)*(detextraterm$sign)
  res=as.numeric(res)

  return(res)

}


####construct the covariance matrix

covStruct=function(theta,cc,dist)
{
  sigma_sq=exp(theta[1])
  alpha=exp(theta[2])

  cov=exp(-dist/alpha)

  res=sigma_sq*cov

  return(res)

}

####construct the covariance matrix

covStruct.nocor=function(theta, n)
{
  sigma_sq=exp(theta)

  return(sigma_sq*diag(rep(1, n)))
}


###get the current estimates of  b using iwls

iwls=function(X, y, D, offset, beta, b, tol, dist)
{

  eta = X%*%beta + b + offset
  mu = exp(eta)

  convEta=10

  while(convEta>tol)
  {
    ##define the working response and weight

    mu=as.vector(mu)
    Y = eta + (y-mu)/mu - offset
    V = diag(1/mu) + D



    solveV=solve(V)


    ###update beta and b
    XVX = crossprod(X, solveV)%*%X
    XVY = crossprod(X, solveV)%*%Y
    DVX = crossprod(D, solveV)%*%X
    DVY = crossprod(D, solveV)%*%Y


    beta.new=solve(XVX, XVY)
    b.new=DVY-DVX%*%beta.new

    eta.new=X%*%beta.new + b.new + offset

    ###covergence
    convEta <- sum((eta.new-eta)^2)/sum(eta.new^2)

    ###updata the eta and mu
    eta = eta.new
    mu = exp(eta)

    #print(convEta)
  }

  return(list(beta=beta.new, b=b.new, eta=eta, mu=mu))

}

###get the current estimates of  b using iwls

iwls_b=function(X, y, D, offset, beta, b, tol, dist)
{

  eta = X%*%beta + b + offset
  mu = exp(eta)
  convEta=10
  solveD=solve(D)



  while(convEta>tol)
  {
    ##define the working response and weight

    mu=as.vector(mu)
    fir_devb = (mu-y) + solveD%*%b
    sec_devb = diag(mu) + solveD

    ###b

    b.new=b-solve(sec_devb)%*%fir_devb
    eta.new=X%*%beta + b.new + offset

    ###covergence
    convEta <- sum((eta.new-eta)^2) /sum(eta^2)

    ###updata the eta and mu
    eta = eta.new


    mu = exp(eta)
    b = b.new

    #print(convEta)
  }

  return(list(beta=as.numeric(beta), b=as.numeric(b), eta=as.numeric(eta), mu=as.numeric(mu)))

}

#function used in adaptive enet

soft.thresholding=function(z, r)
{
   res=0

   if(z>0 &r<abs(z))
   {
     res=z-r
   }

   if(z<0 &r<abs(z))
   {
     res=z+r
   }

  return(res)
}

#function for adaptive elastic net, used in SpatialVS

weighted.ada.enet=function(xmat, y, obs.wts, alpha, lambda, pf1, pf2)
{
  max.iter=100
  eps=1e-4

  wts.sqrt=sqrt(obs.wts)
  xmat.wtd=sweep(xmat, 1, wts.sqrt, "*")
  y.wtd=y*wts.sqrt

  xx.mat=t(xmat.wtd)%*%xmat.wtd
  xx.mat.inv=solve(xx.mat)
  beta0=xx.mat.inv%*%t(xmat.wtd)%*%y.wtd

  beta0=as.vector(beta0)
  beta1=beta0

  pp=length(beta0)

  aa=colSums(xmat.wtd^2)+2*lambda*(1-alpha)*pf2
  aa=as.vector(aa)
  cc=pf1*lambda*alpha

  max.dd=1
  iter=1


  while(iter<max.iter & max.dd>eps)
  {
    for(j in 1:pp)
    {
      rr=y.wtd-xmat.wtd%*%beta1
      dyy=rr+xmat.wtd[,j]*beta1[j]
      bb=sum(dyy*xmat.wtd[,j])
      btmp=soft.thresholding(z=bb, r=cc[j])
      beta1[j]=btmp/aa[j]
    }

    iter=iter+1
    max.dd=max(abs(beta1-beta0))
    beta0=beta1

    #cat(max.dd, iter, "\n")

  }



  res=list(coef=beta1)

  return(res)
}

# calculate initial values for the simulation study

spatialVS_simu_ini=function(simu.data, intercept=T, verbose=T)
{
  y=simu.data$y
  X=simu.data$X
  dist=simu.data$dist
  offset=simu.data$offset

  ID=1:length(y)
  glmmPQL.ini=glmmPQL(y~X-1+offset(offset), random=~1|ID, family=poisson, verbose = F)

  sigma.ini=try(exp(attr(glmmPQL.ini$apVar, "Pars")["reStruct.ID"]), silent=T)

  b.ini=as.numeric(ranef(glmmPQL.ini)[,1])
  beta.ini=fixef(glmmPQL.ini)

  if(class(sigma.ini)=="try-error"|sigma.ini<0.001)
  {
    sigma.ini=0.3
    b.ini=rnorm(length(y), mean=0, sd=sigma.ini)

    D=covStruct(theta=log(c(sigma.ini^2, 4)), cc=0, dist=dist)
    b.ini=as.numeric(iwls_b(X=X, y=y, D=D, offset=offset, beta=beta.ini, b=b.ini, tol=0.001, dist=dist)$b)
  }

  dd.ini=optim(par=log(1.5), fn=cov.fun.fixsigma, method="Brent", lower=-5, upper=10, sigma.ini=sigma.ini, dist=dist, glmmPQL.ini=glmmPQL.ini,beta.ini=beta.ini, b.ini=b.ini, y=y, X=X)$par

  start.ini=list(beta=beta.ini, theta=c(log(sigma.ini^2), dd.ini), b=b.ini)

  tmp=SpatialVS_P_testPQL(data=simu.data, start=start.ini, alpha.vec=0, lambda.vec=0, cc=0, adaptive=F, intercept=intercept, verbose=verbose)

  start=list(beta=as.numeric(start.ini$beta), theta=as.vector(start.ini$theta), b=as.numeric(start.ini$b), lasso.weight=as.numeric(start.ini$beta))


  return(start=start)
}


#Used by SpatialVS in the estimation of covariance paramters under spatial correlation

cov.fun.fixsigma=function(par, sigma.ini, dist, glmmPQL.ini, beta.ini,  b.ini, y, X)
{

  D=covStruct(theta=c(log(sigma.ini^2), par),cc=0, dist=dist)

  mu=exp(glmmPQL.ini$fitted[,"fixed"]+b.ini)

  V=diag(1/mu)+D

  Vinv=solve(V)

  ytrans=b.ini+(y-mu)/mu

  return(as.numeric(determinant(V)$modulus+determinant(t(X)%*%Vinv%*%X)$modulus+t(ytrans-X%*%beta.ini)%*%Vinv%*%(ytrans-X%*%beta.ini)))
}


#intermediate function, performs iterated weighted least squares

iwls_b_new=function(X, y, D, offset, beta, b, tol, dist)    #penalty for zero sum
{
  #browser()

  eta = X%*%beta + b + offset
  mu = exp(eta)
  convEta=10
  solveD=solve(D)

  nn=length(b)
  zmat=matrix(1, ncol=nn, nrow=nn)
  #b=rep(0, nn)

  nIter=0

  while(convEta>tol)
  {
    ##define the working response and weight

    nIter=nIter+1

    mu=as.vector(mu)
    fir_devb = (mu-y) + solveD%*%b+zmat%*%b
    sec_devb = diag(mu) + solveD+ zmat

    ###b
    #browser()
    b.new=b-solve(sec_devb)%*%fir_devb
    eta.new=X%*%beta + b.new + offset

    ###covergence
    convEta <- sum((eta.new-eta)^2) /sum(eta^2)

    #cat("convEta=", convEta, "nIter=", nIter, "mean.b=", mean(b.new), "\n")

    ###updata the eta and mu
    eta = eta.new
    #browser()

    mu = exp(eta)
    b = b.new

    #print(convEta)
  }

  #cat("\n")

  return(list(beta=as.numeric(beta), b=as.numeric(b), eta=as.numeric(eta), mu=as.numeric(mu)))

}


