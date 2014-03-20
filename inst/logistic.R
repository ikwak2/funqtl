## library for the matern covariance function
library(geoR)

###########################################
## functions to generate covariance matrices
###########################################

## creates covariance matrix from the covariance function
## u=time points
## phi=scale
## kappa=smoothness
covMatern <- function(u,phi,kappa)
{
  toeplitz(matern(u,phi,kappa))
}

## creates covariance matrix for equicorrelated variables
## p = dimension of the random vector
## r = common correlation
covEquicor <- function(p,r)
{
  JJ <- matrix(1.0,ncol=p,nrow=p) * r
  diag(JJ) <- 1.0
  return(JJ)
  ##II <- diag(rep(1,p))
  ## we do not check for singularity
  ##(1-r) * II + r * JJ
}


## creates covariance matrix for autocorrelated variables
## u = time points
## r = correlation between two points separated by unit distance
covAutocor <- function(u,r)
{
  toeplitz(exp(log(r)*abs(u-u[1])))
}


#################################################################
## functions to simulate multivariate normal with given covariance
#################################################################

## a potential replacement for next three functions
sim.gaus <- function(n, s2=1., cov.fcn=covAutocor, ...){
  sigma <- s2 * cov.fcn(...)
  ee <- rnormMulti(n, sigma)
  return(ee)
}

rnormMatern <- function(n,u,phi,kappa,s2=1)
  {
    sigma <- s2*covMatern(u,phi,kappa)
    ee <- rnormMulti(n,sigma)
    ee
  }

rnormEquicor <- function(n,p,r,s2=1)
  {
    sigma <- s2*covEquicor(p,r)
    ee <- rnormMulti(n,sigma)
    ee
  }

rnormAutocor <- function(n,u,r,s2=1)
  {
    sigma <- s2*covAutocor(u,r)
    ee <- rnormMulti(n,sigma)
    ee
  }

##########################################
## function to generate multivariate normal
##########################################

## function to simulate a multivariate normal distribution with
## an arbitrary covariance matrix (and zero mean)
rnormMulti <- function(n,sigma)
  {
    ## d <- dim(sigma)
    ## if(d[1]!=d[2])
    ##   stop("Sigma not square.")
    ## else
    
    ## chol() checks for square size in addition to positive-definiteness
    p <- nrow(sigma)
    ee <- rnorm(n*p)
    A <- Matrix::chol(sigma)
    ee <- matrix(ee,nrow=n)%*%A
    
    return(ee)
  }

################################
## logistic growth curve function
################################

## logistic growth curve as a function
logisticFun <- function(beta,tt)
  {
    ## a/(1+b*exp(-c*t))
    ## initial value = a/(1+b)
    ## asymptotic limit = a
    ## rate of growth = c
    beta[1] / (1+ beta[2]*exp(-beta[3]*tt))
  }


##################################################################
## simulates logistic functional data with autocorrelated error
##################################################################
## why is this not using multivariate functions above?
logisticSim1 <- function(beta,tt,n,df=0)
  {
    ## autocorrelation
    rr <- beta[4]
    ## residual variance
    ss <- beta[5]

    m <- length(tt)

    ww <- rr/sqrt(1-rr^2)
    
    mm <- logisticFun(beta[1:3],tt)
    mm <-  kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )
    if(df==0)
      {
        ee <- matrix(rnorm(n*m),nrow=n,ncol=m)
      }
    else
      {
        scaleFac <- sqrt(df/(df-2))
        ee <- matrix(rt(n*m,df=df),nrow=n,ncol=m)/scaleFac
      }
    for( i in 2:m )
      {
        ee[,i] <- (ww*ee[,i-1] + ee[,i])/sqrt(1+ww^2)
      }
    mm + sqrt(ss)*ee
  }

##################################################################
## simulates logistic functional data with Matern error
##################################################################

logisticSim1Mat <- function(beta,tt,n,phi,kappa)
{
  ## residual variance
  ss <- beta[5]

  m <- length(tt)

  mm <- logisticFun(beta[1:3],tt)
  mm <-  kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )


  ee <- rnormMatern(n,tt,phi,kappa)
  
  mm + sqrt(ss)*ee
}

###########################################################################
## Group simulation of logistic functional data with autocorrelated error 
###########################################################################
## abstract group simulation from two ...Sim2 functions
group.sim.logist <- function(beta, tt, grp, single.sim=logisticSim1, ...){
  ngrp <- length(table(grp))# number of groups
  grp <- as.numeric(as.factor(grp))-1

  y <- matrix(nrow=length(grp),ncol=length(tt))
  
  if( nrow(beta)!= ngrp )
    stop("Number of groups does not match the dimension of parameters.")

  for( i in 1:ngrp )
    {
      idx <- which(grp==(i-1))
      y[idx,] <- single.sim(beta[i,],tt,length(idx),...)
    }

  return(y)
}
logisticSim2 <- function(beta,tt,grp,df=0)
{
  ## grp is a vector denoting the group membership
  ngrp <- length(table(grp))# number of groups
  grp <- as.numeric(as.factor(grp))-1

  y <- matrix(nrow=length(grp),ncol=length(tt))
  
  if( nrow(beta)!= ngrp )
    stop("Number of groups does not match the dimension of parameters.")

  for( i in 1:ngrp )
    {
      idx <- which(grp==(i-1))
      y[idx,] <- logisticSim1(beta[i,],tt,length(idx),df=df)
    }

  return(y)
}

###########################################################################
## Group simulation of logistic functional data with Matern error 
###########################################################################

logisticSim2Mat <- function(beta,tt,grp,phi,kappa)
  ## beta :  a matrix of simulation parameters to be passed to simulation functions.
{    
  ngrp <- length(table(grp)) # number of groups
  grp <- as.numeric(as.factor(grp))-1 # why -1?
  y <- matrix(nrow=length(grp),ncol=length(tt))
  
  if( nrow(beta)!= ngrp )
    stop("Number of groups does not match the dimension of parameters.")
    
  for( i in 1:ngrp )
    {
      idx <- which(grp==(i-1))
      y[idx,] <- logisticSim1Mat(beta[i,],tt,length(idx),
                                 phi=phi,kappa=kappa)
    }
  
  return(y)
}

###################################################################
## function to calculate the log likelihood function for functional
## data with logistic mean function and autocorrelated error
###################################################################

logisticLik <- function(beta,y,tt)
{
  rr <- beta[4]
  ss <- beta[5]

  n <- nrow(y)
  m <- ncol(y)

  mm <- logisticFun(beta[1:3],tt)   # mean
  ee <- y - kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )   # residuals

  exponent <- sum((ee[,-m]-rr*ee[,-1])^2)/(1-rr^2) + sum(ee[,m]^2)
  exponent <- (-1./2.)*exponent/ss
  
  loglik <- exponent - m*n*(1/2)*log(2*pi*ss) - n*(m-1)*(1/2)*log(1-rr^2)

  return(loglik)
}

#######################################################################
## function to calculate log likelihood of functional data with
## logistic mean function and autocorrelated error; the parameters have
## been transformed to make the scale better for optimzation function
#######################################################################

logisticLik1 <- function(beta,y,tt)
{
  beta1 <- c(exp(beta[1:3]),tanh(beta[4]),exp(beta[5]))
  return(-logisticLik(beta1,y,tt))
}

## calculation of mle with the 

logisticMLE1 <- function(y,tt,beta0=c(0,0,0,0,0),loglik=FALSE,...)
{
  out <- optim(beta0,logisticLik1,y=y,tt=tt,control=list(fnscale=-1),
               method="BFGS")
  betahat <- c(exp(out$par[1:3]),tanh(out$par[4]),exp(out$par[5]))
  if(loglik)
    {
      loglik <- out$value
      val <- list(betahat,loglik)
    }
  else
    val <- betahat

  val
}



logisticLik2 <- function(beta,y,tt,grp)
{
  grp <- as.numeric(as.factor(grp))-1
  ngrp <- max(grp)+1

  if(length(beta)!=ngrp*3+2)
    stop("Wrong parameter vector length.")
  else
    {
      beta1 <- matrix(nrow=ngrp,ncol=5)
      beta1[,1:3] <- matrix(beta[1:(ngrp*3)],nrow=ngrp,byrow=T)
      beta1[,4:5] <- matrix(rep(beta[(ngrp*3+1):(ngrp*3+2)],ngrp),
                            nrow=ngrp,byrow=T)
      beta1[,c(1,2,3,5)] <- exp(beta1[,c(1,2,3,5)])
      beta1[,4] <- tanh(beta1[,4])
    }
  ## print(beta1)
  beta <- beta1
  loglik <- 0
  for( i in 1:ngrp )
    {
      idx <- which(grp==(i-1))
      loglik <- loglik + logisticLik(beta[i,],y[idx,],tt)
      ## print(logisticLik(beta[i,],y[idx,],tt))
    }
  -loglik
}


logisticMLE2 <- function(y,tt,grp,beta0=NULL,
                         loglik=FALSE,...)
{
  ngrp <- length(unique(grp))

  if(is.null(beta0))
    {
      aa <- max(y)
      if(min(y)<=0)
        bb <- mean(y[,1])
      else
        bb <- (aa-min(y))/min(y)
      cc <- 1
      rr <- 0
      ss <- mean(diag(var(y)))
      beta0 <- c(log(aa),log(bb),log(cc),atanh(rr),log(ss))
      beta0 <- c(rep(beta0[1:3],ngrp),beta0[4:5])
    }

  ## out <- optim(beta0,logisticLik2,y=y,tt=tt,grp=grp,...)
  ## out <- nlm(logisticLik2,beta0,y=y,tt=tt,grp=grp,...)
  out <- nlminb(beta0,logisticLik2,y=y,tt=tt,grp=grp,...)
  ## out$par <- out$estimate
  ## out$convergence <-  as.numeric(out$code>2)
  ## out$value <- out$minimum
  out$value <-  out$objective

  betahat <- matrix(nrow=ngrp,ncol=5)
  betahat[,1:3] <- matrix(out$par[1:(ngrp*3)],ncol=3,byrow=T)
  betahat[,4:5] <- matrix(rep(out$par[(ngrp*3+1):(ngrp*3+2)],ngrp),
                          ncol=2,byrow=T)
  betahat[,c(1,2,3,5)] <- exp(betahat[,c(1,2,3,5)])
  betahat[,4] <- tanh(betahat[,4])
  if(loglik)
    {
      loglik <- out$value
      val <- list(betahat=betahat,loglik=loglik,convergence=out$convergence)
    }
  else
    val <- list(betahat=betahat,convergence=out$convergence)

  val
}

logisticTest2 <- function(y,tt,grp,beta0=NULL,...)
  {
    out1 <- logisticMLE2(y,tt,grp,loglik=TRUE,beta0=beta0,...)

    out0 <- logisticMLE2(y,tt,grp=rep(0,length(grp)),loglik=TRUE,
                         beta0=apply(out1$betahat,2,mean),...)
    
    if((out0$convergence!=0)|(out1$convergence!=0))
      {
        stop("Did not converge.")
        ## out <- list(beta0=NA,beta1=NA,loglik=NA)
      }
    else
      {
        out <- list(beta0=out0$betahat,
                    beta1=out1$betahat,
                    loglik= out0$loglik - out1$loglik)
      }
    return(out)
}


## comparison of methods using data with autoregressive errors
compareSimAR <- function(beta,nsim,tt,psi,grp,df=0)
{    
  ## pval <- matrix(nrow=nsim,ncol=3)
  pval <- matrix(nrow=nsim,ncol=5)    
  n <- length(grp)
  ngrp <- length(unique(grp))
  psidf <- ncol(model.matrix(~psi))*(ngrp-1)
  
  for( i in 1:nsim )
    {
      print(i)
      y <- logisticSim2(beta,tt,grp,df=df)
      out1 <- try(logisticTest2(y,tt,grp,
                                beta0=c(beta[1,1:3],beta[2,1:3],
                                  beta[1,4:5])),TRUE)
      out2 <- quadForm(y,grp,psi,shrink=TRUE)
      out3 <- quadForm(y,grp,psi,shrink=FALSE)

      if(length(grep("loglik",names(out1)))!=0)
        pval[i,1] <- pchisq(2*out1$loglik,3,lower=FALSE)
      else
        pval[i,1] <- NA
      pval[i,2] <- pchisq(out2$quadform,psidf,lower.tail=FALSE)
      pval[i,3] <- pchisq(out3$quadform,psidf,lower.tail=FALSE)        
      pval[i,4] <- photel(out2$quadform,n-ngrp,psidf,lower.tail=FALSE)
      pval[i,5] <- photel(out3$quadform,n-ngrp,psidf,lower.tail=FALSE)        
    }
  return(pval)
}

## comparison of methods using data with Matern errors
compareSimMat <- function(beta,nsim,tt,psi,grp,phi,kappa)
{
  ## pval <- matrix(nrow=nsim,ncol=3)
  pval <- matrix(nrow=nsim,ncol=5)    
  n <- length(grp)
  ngrp <- length(unique(grp))
  psidf <- ncol(model.matrix(~psi))*(ngrp-1)
  
  for( i in 1:nsim )
    {
      print(i)
      y <- logisticSim2Mat(beta,tt,grp,phi=phi,kappa=kappa)
      out1 <- try(logisticTest2(y,tt,grp,
                                beta0=c(beta[1,1:3],beta[2,1:3],
                                  beta[1,4:5])),TRUE)
      out2 <- quadForm(y,grp,psi,shrink=TRUE)        
      out3 <- quadForm(y,grp,psi,shrink=FALSE)
      
      if(length(grep("loglik",names(out1)))!=0)
        pval[i,1] <- pchisq(2*out1$loglik,3,lower=FALSE)
      else
        pval[i,1] <- NA
      pval[i,2] <- pchisq(out2$quadform,psidf,lower.tail=FALSE)
      pval[i,3] <- pchisq(out3$quadform,psidf,lower.tail=FALSE)        
      pval[i,4] <- photel(out2$quadform,n-ngrp,psidf,lower.tail=FALSE)
      pval[i,5] <- photel(out3$quadform,n-ngrp,psidf,lower.tail=FALSE)        
    }
  return(pval)
}

####################################################################

## kolmogorov-smirnov test applied to each column of p-values
ksTest <- function(x)
{
  z <- vector(length=ncol(x))
  for( i in 1:ncol(x) )
    z[i] <- ks.test(x[,i],"punif")$p.value
  return(z)
}

powerCalc <- function(x0,x1,alpha=0.05)
{
  z <- vector(length=col(x0))
  for( i in 1:ncol(x0) )
    z[i] <- mean(x1[,i]<quantile(x0[,i],alpha,na.rm=T),na.rm=T)
  z
}

## function to calculate empirical type-I error for selected colums
sizeCalc <- function(x,size=c(0.1,0.05,0.01),cols=c(1,4,5))
  {
    z <- matrix(nrow=length(size),ncol=length(cols))
    for( j in 1:length(cols) )
      for( i in 1:length(size) )
        {
          z[i,j] <- mean(x[,cols[j]]<size[i],na.rm=T)
        }
    z
  }

###############################
## ROC
###############################

plotROC <- function(x0,x1,cols=1:3,lty=1:length(cols))
  {
    plot(c(0,1),c(0,1),type="n",xaxs="i",yaxs="i",ylim=c(0,1),xlim=c(0,1))
    p <- seq(0,1,by=0.01)
    for( i in 1:length(cols) )
      {
        cutoffs <- quantile(x0[,cols[i]],p,na.rm=T)
        cutoffs[1] <- 0
        cutoffs[length(p)] <- 1
        sens <- c(0,cumsum(hist(x1[,cols[i]],breaks=cutoffs,plot=FALSE)$counts))
        lines(p,sens/max(sens),lty=lty[i])
      }
  }

##############################
## hotelling's t^2 functions
##############################

## we use the fact that T^2(m,p) is (m*p/(m-p+1)) F(p,m-p+1)

photel <- function(x,m,p,lower.tail=TRUE,log.p=FALSE)
  {
    scaleFac <- m*p/(m-p+1)
    pf(x/scaleFac,p,m-p+1,lower.tail=lower.tail,log.p=log.p)
  }

dhotel <- function(x,m,p)
  {
    scaleFac <- m*p/(m-p+1)
    df(x/scaleFac,p,m-p+1)/scaleFac
  }

