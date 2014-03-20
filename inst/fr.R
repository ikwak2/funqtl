#########################################
# Routines for functional regression
#########################################

# y = data matrix
# z = covariate matrix
# phi = smoothing basis

library(qtl)
library(splines)
library(Matrix)
library(corpcor)

varEst <- function(y)
  {
    var(t(y))
  }

bHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
  {
    ty <- t(y)
    if(is.null(phi))
      out1 <- lm(ty~1,weights=weightPhi)
    else
      {
        if(addPhiIntercept)
          out1 <- lm(ty~phi,weights=weightPhi)
        else
          out1 <- lm(ty~phi-1,weights=weightPhi)
      }
    tb <- coef(out1)
    ## print(dim(t(tb)))
    ## print(tb)
    if(is.null(z))
      out2 <- lm(t(tb)~1,weights=weightZ)
    else
      {
        if(addZIntercept)
          out2 <- lm(t(tb)~z,weights=weightZ)
        else
          out2 <- lm(t(tb)~z-1,weights=weightZ)
      }
    b <- coef(out2)
    b
  }

varBHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                    shrink=TRUE)
  {
    # estimate yhat
    yhat <- yHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept)

    if(is.null(z))
      z <- rep(1,length=nrow(y))
    else
      {
        if(addZIntercept)
          z <- model.matrix(~z)
        else
          z <- model.matrix(~z-1)
      }

    # degrees of freedom
    zdf <- ncol(z)
    # sample size
    n <- nrow(y)
    
    # estimate within sample variance matrix, sigma
    if(shrink)
      {
        sigma <- cov.shrink(y-yhat,verb=FALSE)
      }
    else
      {
        sigma <- var(y-yhat)*(n-1)/(n-zdf)
      }
    # project sigma on phi
    pp <- bHat(sigma,phi,phi,addZIntercept=addPhiIntercept,
               addPhiIntercept=addPhiIntercept)

    # calculate (z'z)^(-1)
    zzinv <- solve(t(z)%*%z)
    # take kron prod
    kronecker(zzinv,pp)
  }

quadForm <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                     notIntercept=TRUE,weightPhi=NULL,shrink=TRUE)
  {
    # estimate yhat
    yhat <- yHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept,weightPhi=weightPhi)
    # estimate within sample variance matrix, sigma
    # sigma <- var(y-yhat)
    # sigma <- cov.shrink(y,verb=FALSE)
    n <- nrow(y)
    if(shrink)
      {
        sigma <- cov.shrink(y-yhat,verb=FALSE)
      }
    else
      {
        if(addZIntercept)
          zdf <- ncol(model.matrix(~z))
        else
          zdf <- ncol(model.matrix(~z-1))
        sigma <- var(y-yhat)*(n-1)/(n-zdf)
      }
    
    # sigma <- var(y)
    # calculate bhat
    b <- bHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept,weightPhi=weightPhi)
    # project sigma on phi
    pp <- bHat(sigma,phi,phi,addZIntercept=addPhiIntercept,
               addPhiIntercept=addPhiIntercept,weightPhi=weightPhi,
               weightZ=weightPhi)
    if(is.null(z))
      z <- rep(1,length=nrow(y))
    else
      {
        if(addZIntercept)
          z <- model.matrix(~z)
        else
          z <- model.matrix(~z-1)
      }
    # calculate (z'z)^(-1)
    zzinv <- solve(t(z)%*%z)
    # take kron prod
    vbhat <- kronecker(zzinv,pp)
    if(notIntercept)
      {
        bb <- c(t(b[-1,]))
        idx <- 1:length(b[1,])
        qq <- solve(vbhat[-idx,-idx],bb)
        qq <- sum(bb*qq)
      }
    else
      {
        bb <- c(t(b))
        qq <- solve(vbhat,bb)
        qq <- sum(bb*qq)
      }
    list(quadform=qq,vbhat=vbhat,b=bb)

  }

yHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL)
  {
    ty <- t(y)
    if(is.null(weightPhi))
      weightPhi <- rep(1,ncol(y))
    
    if(is.null(phi))
      {
        phi <- matrix(rep(1,length=ncol(y)),ncol=1)
        # out1 <- lm(ty~phi,weights=weightPhi)
        tyhat <- lm.wfit(x=phi,y=ty,w=weightPhi)$fitted.values
      }
    else
      {
        if(addPhiIntercept)
          # out1 <- lm(ty~phi,weights=weightPhi)
          tyhat <- lm.wfit(x=cbind(rep(1,length=ncol(y)),phi),
                           y=ty,w=weightPhi)$fitted.values
        else
          # out1 <- lm(ty~phi-1,weights=weightPhi)
          tyhat <- lm.wfit(x=phi,y=ty,w=weightPhi)$fitted.values
      }
    # tyhat <- fitted(out1)
      
    if(is.null(z))
      {
        z <- matrix(rep(1,length=nrow(y)),ncol=1)
        # out2 <- lm(t(tyhat)~z)
        yhat <- lm.fit(x=z,y=t(tyhat))$fitted.values
      }
    else
      {
        if(addZIntercept)
          # out2 <- lm(t(tyhat)~z)
          yhat <- lm.fit(x=cbind(1,nrow(y),z),y=t(tyhat))$fitted.values
        else
          # out2 <- lm(t(tyhat)~z-1)
          yhat <- lm.fit(x=z,y=t(tyhat))$fitted.values
        }
    # yhat <- fitted(out2)
    yhat
  }



ySmoothPhi <- function(y,phi,addPhiIntercept=TRUE,weightPhi=NULL)
  {
    ty <- t(y)
    if(is.null(weightPhi))
      weightPhi <- rep(1,ncol(y))
    
    if(is.null(phi))
      {
        phi <- matrix(rep(1,length=ncol(y)),ncol=1)
        # out1 <- lm(ty~phi,weights=weightPhi)
        tyhat <- lm.wfit(x=phi,y=ty,w=weightPhi)$fitted.values
      }
    else
      {
        if(addPhiIntercept)
          # out1 <- lm(ty~phi,weights=weightPhi)
          tyhat <- lm.wfit(x=cbind(rep(1,length=ncol(y)),phi),
                           y=ty,w=weightPhi)$fitted.values
        else
          # out1 <- lm(ty~phi-1,weights=weightPhi)
          tyhat <- lm.wfit(x=phi,y=ty,w=weightPhi)$fitted.values
      }
    t(tyhat)
  }


yHatSmoothPhi <- function(y,z,addZIntercept=TRUE)
  {
    if(is.null(z))
      {
        z <- matrix(rep(1,length=nrow(y)),ncol=1)
        # out2 <- lm(t(tyhat)~z)
        yhat <- lm.fit(x=z,y=y)$fitted.values
      }
    else
      {
        if(addZIntercept)
          # out2 <- lm(t(tyhat)~z)
          yhat <- lm.fit(x=cbind(1,nrow(y),z),y=y)$fitted.values
        else
          # out2 <- lm(t(tyhat)~z-1)
          yhat <- lm.fit(x=z,y=y)$fitted.values
        }
    # yhat <- fitted(out2)
    yhat
  }


devSSSmoothPhi <- function(y,z,addZIntercept=TRUE,weightPhi=NULL)
{
  yh <- yHatSmoothPhi(y,z,addZIntercept=TRUE)
  if(!is.null(weightPhi))
    sum((y-yh)^2%*%diag(weightPhi))
  else
    sum((y-yh)^2)
}


devSS <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                  weightPhi=NULL)
{
  yh <- yHat(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
             weightPhi=weightPhi)
  if(!is.null(weightPhi))
    sum((y-yh)^2%*%diag(weightPhi))
  else
    sum((y-yh)^2)
}


funcFit <- function(y,z,phi,type="ss",calcNull=FALSE,addPhiIntercept=TRUE,
                    addZIntercept=TRUE)
  {
    yhat1 <- yHat(y,z,phi)
    ss1 <- sum((y-yhat1)^2)
    if(calcNull)
      {
        yhat0 <- yHat(y,z=NULL,phi)
        ss0 <- sum((y-yhat0)^2)
      }
    else
      {
        ss0 <- NULL
      }

    if(type=="all")
      {
        b <- bHat(y,z,phi)
        Vbhat <- varBHat(y,z,phi)
        out <- list(b=b,yhat=yhat,Vbhat=Vbhat,ss1=ss1,ss0=ss0)
      }
    else
      {
        out <- list(ss1=ss1,ss0=ss0)
      }
    out
  }

## funcScanone <- function(y,geno,phi)
##   {
##     fit0 <- funcFit(y,geno[,1],phi,calcNull=TRUE)
##     out <- rep(fit0$ss0,ncol(geno))
##     for( i in 1:ncol(geno))
##       {
##         # print(i)
##         fit1 <- funcFit(y,geno[,i],phi)
##         out[i] <- log(fit0$ss0) - log(fit1$ss1)
##       }
##     out
##   }

 funcScanone <- function(y,cr,phi,addcovar=NULL,method="hk",crit="qf",
                         weightPhi=NULL,shrink=TRUE)
  {
    if(method=="hk")
      out <- funcScanone.hk(y,cr,phi,addcovar=addcovar,
                            crit=crit,weightPhi=weightPhi,
                            shrink=shrink)
    else
      if(method=="imp")
        out <- funcScanone.imp(y,cr,phi,addcovar=addcovar,
                               crit=crit,weightPhi=weightPhi)
      else
        error("Unknown method.")
    out
  }

 funcScanone.hk <- function(y,cr,phi,addcovar=NULL,crit="qf",weightPhi=NULL,
                            shrink=TRUE)
   {
     if(!match("prob",names(cr$geno[[1]])))
       {
         warning("First running calc.genoprob.")
         cr <- calc.genoprob(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1],method="hk")
     npseudo <- 0
     nr <- nrow(out)
     gg <- array(dim=c(nrow(y),nr,dim(cr$geno[[1]]$prob)[3]))
     
     for( i in 1:length(cr$geno) )
       {
         # all the genotype probabilities on the i-th chromosome
         ggChr <- cr$geno[[i]]$prob
         npseudoChr <- dim(ggChr)[2]
         # print(dim(ggChr))
         # print(c(npseudo,npseudoChr))
         gg[,(npseudo+1):(npseudo+npseudoChr),] <- ggChr
         npseudo <- npseudo+npseudoChr
       }

     if( crit=="qf" )
       {
         for( j in 1:nr )
           {
             z <- cbind(addcovar,gg[,j,-1])
             qq <- quadForm(y,z,phi,weightPhi=weightPhi,
                            shrink=shrink)
             out[j,3] <- qq$quadform/2/log(10)
           }
       }

     if( crit=="ss" )
       {
         z <- addcovar
         ySmooth <- ySmoothPhi(y,phi,weightPhi=weightPhi)
         ss0 <- devSSSmoothPhi(ySmooth,addcovar,weightPhi=weightPhi)
         # ss0 <- devSS(y,addcovar,phi,weightPhi=weightPhi)
         for( j in 1:nr )
           {
             z <- cbind(addcovar,gg[,j,-1])
            ss  <- devSSSmoothPhi(ySmooth,z,weightPhi=weightPhi)
            out[j,3] <- -log(ss/ss0)
           }
       }
     out
   }

 funcScanone.imp <- function(y,cr,phi)
   {
     if(!match("draws",names(cr$geno[[1]])))
       {
         warning("First running sim.geno.")
         cr <- sim.geno(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1],method="imp")
     
     npseudo <- 0
     for( i in 1:length(cr$geno) )
       {
         # all the imputations on the i-th chromosome
         gg <- cr$geno[[i]]$draws
         npseudoChr <- dim(gg)[2]
         for( j in 1:npseudoChr )
           {
             qq <- quadForm(y,gg[,j],phi)
             out[npseudo+j,3] <- qq$quadform/2/log(10)
           }
         npseudo <- npseudo+npseudoChr
       }
     out
   }


funcScanonePerm <- function(y,cr,phi,nperm,addcovar=NULL,method="hk",crit="qf",
                            weightPhi=NULL,shrink=TRUE)
  {
    # if genotype probabilities have not been calculated, calculate them
     if(!match("prob",names(cr$geno[[1]])))
       {
         warning("First running calc.genoprob.")
         cr <- calc.genoprob(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1])
     out[,3:(2+nperm)] <- matrix(rep(0,nperm*nrow(out)),ncol=nperm)

     # number of individuals
     n <- nrow(y)
     m <- ncol(y)
     
     # ready genotype matrix
     npseudo <- 0
     # number pseudomarkers
     nr <- nrow(out)
     gg <- array(dim=c(nrow(y),nr,dim(cr$geno[[1]]$prob)[3]))
     
     for( i in 1:length(cr$geno) )
       {
         # all the genotype probabilities on the i-th chromosome
         ggChr <- cr$geno[[i]]$prob
         npseudoChr <- dim(ggChr)[2]
         # print(dim(ggChr))
         # print(c(npseudo,npseudoChr))
         gg[,(npseudo+1):(npseudo+npseudoChr),] <- ggChr
         npseudo <- npseudo+npseudoChr
       }


     if( crit=="ss" )
       {
         ss0 <- devSS(y,addcovar,phi,weightPhi=weightPhi)
       }
     
     # perform permutations
     for( k in 1:nperm )
       {
         print(k)
         npseudo <- 0
         permidx <- sample(n)
         if( crit=="qf" )
           {
             for( j in 1:nr )
               {
                 qq <- quadForm(y[permidx,],gg[,j,-1],phi,weightPhi=weightPhi,
                                shrink=shrink)
                 out[j,2+k] <- qq$quadform/2/log(10)
           }
           }
         
         if( crit=="ss" )
           {
             for( j in 1:nr )
               {
                 if(!is.null(addcovar))
                   {
                     z <- cbind(addcovar[permidx,],gg[,j,-1])
                   }
                 else
                   {
                     z <- gg[,j,-1]
                   }
                 ss  <- devSS(y[permidx,],z,phi,weightPhi=weightPhi)
                 out[j,2+k] <- -log(ss/ss0)
               }
           }
       }
     # to conform to R/qtl scanone.perm
     out <- apply(out[,-(1:2)],2,max)
     out <- matrix(out,ncol=1)
     class(out) <- c("scanoneperm", "matrix")
     out
  }

funcScanonePermCluster <- function(y,cr,phi,nperm,addcovar=NULL,
                                   method="hk",crit="qf",weightPhi=NULL,
                                   shrink=TRUE,ncluster)
  {
    cl <- makeCluster(ncluster)
    clusterStopped <- FALSE
    on.exit(if (!clusterStopped) stopCluster(cl))
    clusterSetupRNG(cl)
    clusterEvalQ(cl, require(qtl, quietly = TRUE))
    clusterEvalQ(cl, source("~/HG/functionalQTL/programs/R/fr.R"))
    nperm <- ceiling(nperm/ncluster)
    # if (missing(chr))
    #  chr <- names(cross$geno)
    operm <- clusterCall(cl, funcScanonePerm,y,cr,phi,nperm,addcovar,
                         method,crit,weightPhi,shrink)
    stopCluster(cl)
    clusterStopped <- TRUE
    for (j in 2:length(operm)) operm[[1]] <- c(operm[[1]],
                                               operm[[j]])
    return(operm[[1]])
   
  }

###################################################################
# pairscan routines
###################################################################

 funcScantwo <- function(y,cr,phi,addcovar=NULL,method="hk",crit="qf",
                         weightPhi=NULL,shrink=TRUE)
  {
    if(method=="hk")
      out <- funcScantwo.hk(y,cr,phi,addcovar=addcovar,
                            crit=crit,weightPhi=weightPhi,
                            shrink=shrink)
    else
      if(method=="imp")
        out <- funcScantwo.imp(y,cr,phi,addcovar=addcovar,
                               crit=crit,weightPhi=weightPhi)
      else
        error("Unknown method.")
    out
  }

 funcScantwo.hk <- function(y,cr,phi,addcovar=NULL,crit="qf",weightPhi=NULL,
                            shrink=TRUE)
   {
     if(!match("prob",names(cr$geno[[1]])))
       {
         warning("First running calc.genoprob.")
         cr <- calc.genoprob(cr)
       }
     # get the array ready
     out1 <- scanone(cr,pheno.col=y[,1],method="hk")
     npseudo <- nrow(out1)
     # make array of genotype probabilities
     gg <- array(dim=c(nrow(y),npseudo,dim(cr$geno[[1]]$prob)[3]))

     out <- matrix(nrow=npseudo,ncol=npseudo)
     
     # make genotype probability matrix
     npseudoCount <- 0
     for( i in 1:length(cr$geno) )
       {
         # all the genotype probabilities on the i-th chromosome
         ggChr <- cr$geno[[i]]$prob
         npseudoChr <- dim(ggChr)[2]
         # print(dim(ggChr))
         # print(c(npseudo,npseudoChr))
         gg[,(npseudoCount+1):(npseudoCount+npseudoChr),] <- ggChr
         npseudoCount <- npseudoCount+npseudoChr
       }

     if( crit=="qf" )
       {
         for( i in 1:(npseudo-1) )
           {
             # first marker with additive covariates
             z1 <- gg[,i,-1]
             zMain <- cbind(addcovar,z1)
             qqMain <- quadForm(y,zMain,phi,weightPhi=weightPhi,
                                shrink=shrink)
             out[i,i] <- qqMain$quadform/2/log(10)
             for( j in (i+1):npseudo)
               {
                 # second marker sans covariate
                 z2 <- gg[,j,-1]
                 zAdd <- cbind(addcovar,z1,z2)
                 qqAdd <- quadForm(y,zAdd,phi,weightPhi=weightPhi,
                                shrink=shrink)
                 out[i,j] <- qqAdd$quadform/2/log(10)
                 zFull <- cbind(zAdd,kronecker(z1,z2))
                 qqFull <- quadForm(y,zFull,phi,weightPhi=weightPhi,
                                shrink=shrink)
                 out[j,i] <- qqFull$quadform/2/log(10)
               }
           }
       }

     if( crit=="ss" )
       {
         # get the smoothed trajectory
         ySmooth <- ySmoothPhi(y,phi,weightPhi=weightPhi)
         # baseline sum of squares
         ss0 <- devSSSmoothPhi(ySmooth,addcovar,weightPhi=weightPhi)

         print(dim(gg))
         for( i in 1:(npseudo-1) )
           {
             print(i)
             # first marker with additive covariates
             z1 <- as.matrix(gg[,i,-1])
             zMain <- cbind(addcovar,z1)
             ssMain  <- devSSSmoothPhi(ySmooth,zMain,weightPhi=weightPhi)
             out[i,i] <- -log(ssMain/ss0)

             # loop over pairs
             for( j in (i+1):(npseudo))
               {
                 # second marker sans covariate
                 z2 <- as.matrix(gg[,j,-1])
                 zAdd <- cbind(addcovar,z1,z2)
                 ssAdd <- devSSSmoothPhi(ySmooth,zAdd,weightPhi=weightPhi)
                 out[i,j] <- -log(ssAdd/ss0)
                 zFull <- cbind(zAdd,model.matrix(~z1:z2-1))
                 ssFull <- devSSSmoothPhi(ySmooth,zFull,weightPhi=weightPhi)
                 out[j,i] <- -log(ssFull/ss0)
               }
           }
       }
#         z <- addcovar
#         ySmooth <- ySmoothPhi(y,phi,weightPhi=weightPhi)
#         ss0 <- devSSSmoothPhi(ySmooth,addcovar,weightPhi=weightPhi)
#         # ss0 <- devSS(y,addcovar,phi,weightPhi=weightPhi)
#         for( j in 1:nr )
#           {
#             z <- cbind(addcovar,gg[,j,-1])
#            ss  <- devSSSmoothPhi(ySmooth,z,weightPhi=weightPhi)
#            out[j,3] <- -log(ss/ss0)
#           }

     map <- data.frame(chr=out1[,1],pos=out1[,2],eq.spacing=1,xchr=FALSE)
     val <- list(lod=out,map=map,scanoneX=NULL)
     class(val) <- "scantwo"
     val
   }

###########################################################
# routine for cross validated integrated sequared error
###########################################################

cvSS <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10)
  {
    # divide the data into nfolds
    n <- nrow(y)
    idx <- sample(n)
    folds <- ((0:(n-1))%%nfolds) + 1
    
    if(addPhiIntercept)
      phi1 <- model.matrix(~phi)
    if(addZIntercept)
      {
        if(is.null(z))
          {
            z1 <- model.matrix(~rep(1,n)-1)
          }
        else
          {
            z1 <- model.matrix(~z)
          }
      }

    ss <- vector(length=nfolds)
    z1 <- as.matrix(z1)
    
    for( i in 1:nfolds )
      {
        # select the i-the fold
        fidx <- idx[folds==i]
        bb <- bHat(y[-fidx,],z1[-fidx,],phi=phi,addPhiIntercept=addPhiIntercept,
                  addZIntercept=FALSE,weightPhi=weightPhi,
                  weightZ=weightZ[-fidx])

        yh <- z1[fidx,]%*%bb%*%t(phi1)

        if(!is.null(weightPhi))
          ss[i] <- sum((y[fidx,]-yh)^2%*%diag(weightPhi))
        else
          ss[i] <- sum((y[fidx,]-yh)^2)
      }
    ss
  }


bsCV <- function(y,z,x,df=3:20,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10,degree=3,
                 intercept=FALSE)
{
  cv <- matrix(nrow=length(df),ncol=2)
  for( i in 1:nrow(cv) )
    {
      phi <- bs(x,df=df[i],degree=degree,intercept=intercept)
      out <- cvSS(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10)
      cv[i,1] <- mean(out)
      cv[i,2] <- sd(out)
    }
  cv
}


yPred <- function(y,z,phi,znew,phinew,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
  {
    bb <- bHat(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
    if(addZIntercept)
      znew <- model.matrix(~znew)
    else
      znew <- model.matrix(~znew-1)

    if(addPhiIntercept)
      phinew <- model.matrix(~phinew)
    else
      phinew <- model.matrix(~phinew-1)

    print(dim(bb))
    
    yp <- znew %*% bb %*% t(phinew)

    yp
  }
  
