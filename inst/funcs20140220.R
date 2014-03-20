library(qtl)
library(fda)

getY <- function(cross, n.max=5, criteria=.9, nn = 0) {

    Y = cross$pheno
    udv <- svd(Y)
    vec <- udv$d^2
    vec <- vec / sum(vec)

    ss = 0

    if( nn == 0) {
        for(j in 1:n.max) {
            ss = ss + vec[j]
            if ( ss > criteria ) {
                nn = j
                break
            }
        }
    }
    pc <- udv$u %*% diag(udv$d)
    pc[,1:nn]
}




getY2 <- function(cross, n.max=4, criteria=.9, nn = 0, basis) {
    Y = t( cross$pheno )
    m = nrow(Y)
#    splinebasis.y <- create.bspline.basis(c(0,m),4,1)
    if(missing(basis)) {
        splinebasis.y <- create.bspline.basis(c(0,m), nbasis = m/2)
    } else {
        splinebasis.y <- basis
    }
#    splinebasis.y <- create.bspline.basis(c(0,m),round(m/2),4)

    time <- 0:(m-1) + 0.5
    mat <- eval.basis(time, splinebasis.y)

    coef.y <- solve(crossprod(mat), crossprod(mat, Y))

    yfd = fd(coef.y, splinebasis.y, list("time","indv","value") )

    if( nn == 0) {
        y.pcalist3 = pca.fd(yfd, n.max)

        ss = 0
        for(j in 1:n.max) {
            ss = ss + y.pcalist3$varprop[j]
            if ( ss > criteria ) {
                nn = j
                break
            }
        }
    }


    y.pcalist3 = pca.fd(yfd, nn)
    eigfc3 <- y.pcalist3$harmonics
    mat3 <- eval.fd(time, eigfc3)
    nY3 <- t(solve(crossprod(mat3), crossprod(mat3, Y) ) )
}



scanoneM2 <- function(cross, Y) {

    n.ind <- nind(cross)
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    resdf <- n.ind - 2
    testdf <- 1
    nt <- 3

    out <- NULL;
    for(i in 1:n.chr) {
        LOD = NULL;
        for(j in 1:n.mar[i]) {
            X <- cross$geno[[i]]$prob[,j,]

            lm1 <- lm(Y ~ X -1)
            B <- lm1$coef
            RSSp <-  t(Y)%*%Y  - t(B) %*% t(X) %*% Y

            X2 <- rep(1,n.ind)
            B2 <- qr.solve(t(X2)%*%X2, t(X2)%*%Y)

            RSSr <- t(Y)%*%Y - t(B2) %*% t(X2) %*% Y
            LOD <- c(LOD, -( resdf - 1/2*(nt-testdf +1) ) * log( det(RSSp)/det(RSSr), 10) )
        }
        out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
    }

    outt <- data.frame( chr = out[,1], pos = out[,2], lod = out[,3])

    class(outt) <- c("scanone", "data.frame")
    outt
}


scanoneM <- function(cross, Y, tol=1e-7) {

    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
#    p = ncol(Y)
    if(is.vector(Y)) { p = 1} else {p = ncol(Y)}

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    L0 <- determinant(Sigma)$modulus


    out <- NULL;
    for(i in 1:n.chr) {
        LOD = NULL;
        for(j in 1:n.mar[i]) {
            X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,])
            E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
            Sigma <- crossprod(E)
            L1 <- determinant(Sigma)$modulus
            LOD <- c(LOD, n.ind/2*log(L0/L1, 10) )
        }
        out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
    }

    outt <- data.frame( chr = out[,1], pos = out[,2], lod = out[,3])

    class(outt) <- c("scanone", "data.frame")
    outt
}


scanoneM3 <- function(cross, Y, tol=1e-7, n.perm=0) {

    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    if(is.vector(Y)) { p = 1} else {p = ncol(Y)}


    if(n.perm == 0) {


        E <- matrix(NA, n.ind, p)
        X <- cbind(rep(1,n.ind))
        E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
        Sigma <- crossprod(E)
        L0 <- determinant(Sigma)$modulus


        out <- NULL;
        for(i in 1:n.chr) {
            LOD = NULL;
            for(j in 1:n.mar[i]) {
                X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
                Sigma <- crossprod(E)
                L1 <- determinant(Sigma)$modulus
#                LOD <- c(LOD, n.ind/2*log(L0/L1,10 ) )
                LOD <- c(LOD, n.ind/2*log10(exp(1))*(L0 - L1) )
            }
            out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
        }

        outt <- data.frame( chr = out[,1], pos = out[,2], lod = out[,3])

        class(outt) <- c("scanone", "data.frame")
        return(outt)
    } else {

        lods = NULL;
        for( rep in 1:n.perm) {
            o <- sample(n.ind)
            if(is.vector(Y)) { nY <- Y[o] } else { nY <- Y[o,] }


            E <- matrix(NA, n.ind, p)
            X <- cbind(rep(1,n.ind))
            E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
            Sigma <- crossprod(E)
            L0 <- determinant(Sigma)$modulus


            out <- NULL;
            for(i in 1:n.chr) {
                LOD = NULL;
                for(j in 1:n.mar[i]) {
                    X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                    E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
                    Sigma <- crossprod(E)
                    L1 <- determinant(Sigma)$modulus
#                    LOD <- c(LOD, n.ind/2*log(L0/L1,10) )
                    LOD <- c(LOD, n.ind/2*log10(exp(1))*(L0 - L1) )

                }
                out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
            }

            lods <- c(lods, max(out[,3]) )
        }
        return(lods)
    }
}



scanoneM <- function(cross, Y, tol=1e-7, n.perm=0, method=c("hk","f", "sl", "ml"), pheno.cols ) {

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)
    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if(n.perm == 0) {

        if (method == "sl" || method == "ml") {
            if (!("prob" %in% names(cross$geno[[1]]))) {
                warning("First running calc.genoprob.")
                cross <- calc.genoprob(cross)
            }
            temp <- cross
            temp$pheno[,1:p] <- Y

            outs <- scanoneF(temp, pheno.cols=1:p)

            return(outs)
        } else {

            E <- matrix(NA, n.ind, p)
            X <- cbind(rep(1,n.ind))
            E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
            Sigma <- crossprod(E)

            if( method == "hk") {
                L0 <- determinant(Sigma)$modulus
            } else {
                L0 <- sum(diag(Sigma))
            }

            out <- NULL;
            for(i in 1:n.chr) {
                LOD = NULL;
                for(j in 1:n.mar[i]) {
                    X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
                    Sigma <- crossprod(E)

                    if( method == "hk") {
                        L1 <- determinant(Sigma)$modulus
                    } else {
                        L1 <- sum(diag(Sigma))
                    }

#                LOD <- c(LOD, n.ind/2*log(L0/L1,10 ) )
                    LOD <- c(LOD, n.ind/2*log10(exp(1))*(L0 - L1) )
                }
                out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
            }

            outt <- data.frame( chr = out[,1], pos = out[,2], lod = out[,3])

            class(outt) <- c("scanone", "data.frame")
            return(outt)
        }
    } else {



        if (method == "sl" || method == "ml") {
            temp <- cross
            temp$pheno[,1:p] <- Y

            outs <- scanoneF(temp, pheno.cols=1:p, n.perm = n.perm)

            return(outs)
        } else {

            lods = NULL;
            for( rep in 1:n.perm) {
                o <- sample(n.ind)
                if(is.vector(Y)) { nY <- Y[o] } else { nY <- Y[o,] }


                E <- matrix(NA, n.ind, p)
                X <- cbind(rep(1,n.ind))
                E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
                Sigma <- crossprod(E)

                if( method == "hk") {
                    L0 <- determinant(Sigma)$modulus
                } else {
                    L0 <- sum(diag(Sigma))
                }

                out <- NULL;
                for(i in 1:n.chr) {
                    LOD = NULL;
                    for(j in 1:n.mar[i]) {
                        X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                        E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
                        Sigma <- crossprod(E)

                        if( method == "hk") {
                            L1 <- determinant(Sigma)$modulus
                        } else {
                            L1 <- sum(diag(Sigma))
                        }
   #                    LOD <- c(LOD, n.ind/2*log(L0/L1,10) )
                        LOD <- c(LOD, n.ind/2*log10(exp(1))*(L0 - L1) )

                    }
                    out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
                }

                lods <- c(lods, max(out[,3]) )
            }
            return(lods)
        }
    }
}




scanoneL <- function(cross, pheno.cols, tol=1e-7, n.perm=0) {

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    n.phe <- nphe(cross)
    Y <- as.matrix(cross$pheno)

    if(n.perm == 0) {


        E <- matrix(NA, n.ind, n.phe)
        X <- cbind(rep(1,n.ind))
        E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
        Sigma <- crossprod(E)
        LMSSE0 <- sum(diag(Sigma))  # sum( diag( t(E) %*% E ) )

        out <- NULL;
        for(i in 1:n.chr) {
            LOD = NULL;
            for(j in 1:n.mar[i]) {
                X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,])
                E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
                Sigma <- crossprod(E)
                LMSSE1 <- sum(diag(Sigma))
#                LOD <- c(LOD, n.ind/2*log(LMSSE0/LMSSE1, 10) )
                LOD <- c(LOD, n.ind/2*log10(exp(1))*(LMSSE0 - LMSSE1) )

            }
            out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
        }

        outt <- data.frame( chr = out[,1], pos = out[,2], lod = out[,3])

        class(outt) <- c("scanone", "data.frame")
        outt
    } else {


        lods = NULL;
        for( rep in 1:n.perm) {
            o <- sample(n.ind)
            if(is.vector(Y)) { nY <- Y[o] } else { nY <- Y[o,] }


            E <- matrix(NA, n.ind, n.phe)
            X <- cbind(rep(1,n.ind))
            E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
            Sigma <- crossprod(E)
            L0 <- sum(diag(Sigma))


            out <- NULL;
            for(i in 1:n.chr) {
                LOD = NULL;
                for(j in 1:n.mar[i]) {
                    X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                    E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
                    Sigma <- crossprod(E)
                    L1 <- sum(diag(Sigma))
#                    LOD <- c(LOD, n.ind/2*log(L0/L1,10) )
                    LOD <- c(LOD, n.ind/2*log10(exp(1))*(L0 - L1) )

                }
                out <- rbind(out, cbind(rep(i,n.mar[i]), cross$geno[[i]]$map, LOD) )
            }

            lods <- c(lods, max(out[,3]) )
        }
        return(lods)
    }
}














getthresholdM <- function(cross, Y, tol=1e-7, n.perm=1000, method = c("hk", "f"), pheno.cols) {

    method <- match.arg(method)

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    lods = NULL;
    for( rep in 1:n.perm) {
        o <- sample(n.ind)
        if(is.vector(Y)) { nY <- Y[o] } else { nY <- Y[o,] }


        E <- matrix(NA, n.ind, p)
        X <- cbind(rep(1,n.ind))
        E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
        Sigma <- crossprod(E)

        if( method == "hk") {
            L0 <- determinant(Sigma)$modulus
        } else {
            L0 <- sum(diag(Sigma))
        }

        LODf = NULL;
        LODa = NULL;
        LOD1 = NULL;
        for(i in 1:n.chr) {
            for(j in 1:n.mar[i]) {

                X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                E <- .Call(stats:::C_Cdqrls, X, nY, tol)$residuals
                Sigmaa <- crossprod(E)


                if( method == "hk") {
                    L1 <- determinant(Sigma)$modulus
                } else {
                    L1 <- sum(diag(Sigma))
                }

                LOD1 <- c(LOD1, n.ind/2*log(L0/L1,10) )

                for(ii in 1:n.chr) {
                    for(jj in 1:n.mar[ii]) {
                        if( (ii == i && j > jj ) && i >= ii ) {
                            Xa <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1], cross$geno[[ii]]$prob[,jj,1])
                            Xf <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1], cross$geno[[ii]]$prob[,jj,1], cross$geno[[i]]$prob[,j,1] * cross$geno[[ii]]$prob[,jj,1])
                            Ea <- .Call(stats:::C_Cdqrls, Xa, nY, tol)$residuals
                            Ef <- .Call(stats:::C_Cdqrls, Xf, nY, tol)$residuals
                            Sigmaa <- crossprod(Ea)
                            Sigmaf <- crossprod(Ef)


                            if( method == "hk") {
                                L1a <- determinant(Sigmaa)$modulus
                                L1f <- determinant(Sigmaf)$modulus
                            } else {
                                L1a <- sum(diag(Sigma))
                                L1f <- sum(diag(Sigma))
                            }

                            LODa <- c(LODa, n.ind/2*log10(exp(1))*(L0 - L1a) )
                            LODf <- c(LODf, n.ind/2*log10(exp(1))*(L0 - L1f) )

                        }
                    }
                }
            }
        }


        lods <- rbind(lods, c(max(LOD1), max(LODa), max(LODf)) )
    }
    colnames(lods) = c("LOD1", "LODa", "LODf")
    return(lods)
}







getlodM2 <- function(cross, Y, formula, qtl, tol=1e-7) {

    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    n.ind <- qtl$n.ind
    n.qtl <- qtl$n.qtl
    n.gen <- qtl$n.gen
    mark <- find.marker(cross, qtl$chr, qtl$pos)
    names(mark) <- qtl$altname
    nx <- length(mark)


    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        formula <- as.formula(formula)
    }
                                        #    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    if(is.vector(Y)) { p = 1} else {p = ncol(Y)}

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    L0 <- determinant(Sigma)$modulus


    tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
    terms <- strsplit(tempform, " *\\+ *")[[1]]

    X <- NULL
    X <- cbind(rep(1,n.ind))

    for (j in 1:length(terms)) {
        Xadd <- NULL;
        if (length(grep(":", terms[j])) > 0) {
            temp <- strsplit(terms[j], " *: *")[[1]]

            wh1 <- which( strsplit(terms[j], " *: *")[[1]][1] == names(mark))
            wh2 <- which( strsplit(terms[j], " *: *")[[1]][2] == names(mark))
            qtt1 <- mark[ wh1 ]
            qtt2 <- mark[ wh2 ]
            i1 <- as.numeric(qtl$chr[ wh1 ])
            i2 <- as.numeric(qtl$chr[ wh2 ])

            Xadd <- cross$geno[[i1]]$prob[ ,qtt1, 1 ] *
                cross$geno[[i2]]$prob[ ,qtt2, 1 ]

        } else {
            wh <- which( terms[j] == names(mark))
            qtt <- mark[ wh ]
            i = as.numeric(qtl$chr[ wh ])

            Xadd <- cross$geno[[i]]$prob[ ,qtt, 1 ]
        }

        X <- cbind(X, Xadd)
    }
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    L1 <- determinant(Sigma)$modulus
    LOD = n.ind/2*log10(exp(1))*(L0 - L1)

}



getlodM <- function(cross, Y, formula, qtl, tol=1e-7, method = c("hk","f"), pheno.cols) {

    method <- match.arg(method)

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }


    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    n.ind <- qtl$n.ind
    n.qtl <- qtl$n.qtl
    n.gen <- qtl$n.gen
    mark <- find.marker(cross, qtl$chr, qtl$pos)
    names(mark) <- qtl$altname
    nx <- length(mark)


    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        formula <- as.formula(formula)
    }
                                        #    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)



    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)


    if( method == "hk") {
        L0 <- determinant(Sigma)$modulus
    } else {
        L0 <- sum(diag(Sigma))
    }



    tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
    terms <- strsplit(tempform, " *\\+ *")[[1]]

    X <- NULL
    X <- cbind(rep(1,n.ind))

# j = 1

    for (j in 1:length(terms)) {
        Xadd <- NULL;
        if (length(grep(":", terms[j])) > 0) {

            if( strsplit(terms[j], " *: *")[[1]][1] %in% names(mark) && strsplit(terms[j], " *: *")[[1]][2] %in% names(mark) ) {
                temp <- strsplit(terms[j], " *: *")[[1]]

                wh1 <- which( strsplit(terms[j], " *: *")[[1]][1] == names(mark))
                wh2 <- which( strsplit(terms[j], " *: *")[[1]][2] == names(mark))
                qtt1 <- mark[ wh1 ]
                qtt2 <- mark[ wh2 ]
                i1 <- (qtl$chr[ wh1 ])
                i2 <- (qtl$chr[ wh2 ])

                Xadd <- cross$geno[[i1]]$prob[ ,qtt1, 1 ] *
                    cross$geno[[i2]]$prob[ ,qtt2, 1 ]
            }
        } else {
            if( terms[j] %in% names(mark) ) {
                wh <- which( terms[j] == names(mark))
                qtt <- mark[ wh ]
                i = (qtl$chr[ wh ])

                Xadd <- cross$geno[[i]]$prob[ ,qtt, 1 ]
            }
        }

        X <- cbind(X, Xadd)
    }
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    if( method == "hk") {
        L1 <- determinant(Sigma)$modulus
    } else {
        L1 <- sum(diag(Sigma))
    }
    LOD = n.ind/2*log10(exp(1))*(L0 - L1)
}



#out <- fitqtlengineM(ril, Y, formula, qtl)


fitqtlengineM <- function(cross, Y, formula, qtl, tol=1e-7, method=c("hk","f"), pheno.cols) {


    method <- match.arg(method)

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }



    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    n.ind <- qtl$n.ind
    n.qtl <- qtl$n.qtl
    n.gen <- qtl$n.gen
    mark <- find.marker(cross, qtl$chr, qtl$pos)
    names(mark) <- qtl$altname
    nx <- length(mark)


    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        formula <- as.formula(formula)
    }
                                        #    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)

    if( method == "hk") {
        L0 <- determinant(Sigma)$modulus
    } else {
        L0 <- sum(diag(Sigma))
    }

    tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
    terms <- strsplit(tempform, " *\\+ *")[[1]]

    X <- NULL
    X <- cbind(rep(1,n.ind))

    for (j in 1:length(terms)) {
        Xadd <- NULL;
        if (length(grep(":", terms[j])) > 0) {
            if( strsplit(terms[j], " *: *")[[1]][1] %in% names(mark) && strsplit(terms[j], " *: *")[[1]][2] %in% names(mark) ) {

                temp <- strsplit(terms[j], " *: *")[[1]]

                wh1 <- which( strsplit(terms[j], " *: *")[[1]][1] == names(mark))
                wh2 <- which( strsplit(terms[j], " *: *")[[1]][2] == names(mark))
                qtt1 <- mark[ wh1 ]
                qtt2 <- mark[ wh2 ]
                i1 <- (qtl$chr[ wh1 ])
                i2 <- (qtl$chr[ wh2 ])

                Xadd <- cross$geno[[i1]]$prob[ ,qtt1, 1 ] *
                    cross$geno[[i2]]$prob[ ,qtt2, 1 ]
            }
        } else {

            if( terms[j] %in% names(mark) ) {
                wh <- which( terms[j] == names(mark))
                qtt <- mark[ wh ]
                i = (qtl$chr[ wh ])

                Xadd <- cross$geno[[i]]$prob[ ,qtt, 1 ]
            }
        }

        X <- cbind(X, Xadd)
    }
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)

    if( method == "hk") {
        L1 <- determinant(Sigma)$modulus
    } else {
        L1 <- sum(diag(Sigma))
    }
    LOD = n.ind/2*log10(exp(1))*(L0 - L1)

    dfM = p * length(terms)
    dfE = p * n.ind - dfM - p
    dfT = p * ( n.ind - 1 )

    SSE = L1
    SST = determinant(crossprod(Y))$modulus
    SSM = SST - SSE


    MSM = SSM / 4
    MSE = SSE / dfE

    result.full <- data.frame( df =c(dfM, dfE,dfT), SS = c(SSM,SSE,SST), MS = c(MSM, MSE,NA), LOD =c( LOD, NA,NA) )
    rownames(result.full) <- c("Model","Error","Total")

    return(list( result.full = result.full, lod = LOD) )

}







#qtl <- makeqtl(ril, chr = c(1, 3, 4),
#               pos = c(64, 17, 40.3), what = "prob")

#formula <- as.formula("y ~ Q1 + Q2 + Q3 + Q1:Q3")

#out <- fitqtlengineM(cross, Y, formula = formula, qtl = qtl)
#out <- fitqtlM(cross, Y, formula, qtl)
#ouout <- fitqtl(cross, pheno.col=1, formula = formula, qtl = qtl, method = "hk")

fitqtlM <- function(cross, Y, formula, qtl, tol=1e-7, method=c("hk","f"), pheno.cols) {

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)


    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    n.ind <- qtl$n.ind
    n.qtl <- qtl$n.qtl
    n.gen <- qtl$n.gen
    mark <- find.marker(cross, qtl$chr, qtl$pos)
    names(mark) <- qtl$altname
    nx <- length(mark)


    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        formula <- as.formula(formula)
    }
                                        #    n.ind <- nind(cross) # n
    n.phe <- nphe(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    if( method == "hk") {
        L0 <- determinant(Sigma)$modulus
    } else {
        L0 <- sum(diag(Sigma))
    }


    tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
    terms <- strsplit(tempform, " *\\+ *")[[1]]

    X <- NULL
    X <- cbind(rep(1,n.ind))

    for (j in 1:length(terms)) {
        Xadd <- NULL;
        if (length(grep(":", terms[j])) > 0) {
            if( strsplit(terms[j], " *: *")[[1]][1] %in% names(mark) && strsplit(terms[j], " *: *")[[1]][2] %in% names(mark) ) {

                temp <- strsplit(terms[j], " *: *")[[1]]

                wh1 <- which( strsplit(terms[j], " *: *")[[1]][1] == names(mark))
                wh2 <- which( strsplit(terms[j], " *: *")[[1]][2] == names(mark))
                qtt1 <- mark[ wh1 ]
                qtt2 <- mark[ wh2 ]
                i1 <- (qtl$chr[ wh1 ])
                i2 <- (qtl$chr[ wh2 ])

                Xadd <- cross$geno[[i1]]$prob[ ,qtt1, 1 ] *
                    cross$geno[[i2]]$prob[ ,qtt2, 1 ]
            }
        } else {

            if( terms[j] %in% names(mark) ) {
                wh <- which( terms[j] == names(mark))
                qtt <- mark[ wh ]
                i = (qtl$chr[ wh ])

                Xadd <- cross$geno[[i]]$prob[ ,qtt, 1 ]
            }
        }

        X <- cbind(X, Xadd)
    }
    E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
    Sigma <- crossprod(E)
    if( method == "hk") {
        L1 <- determinant(Sigma)$modulus
    } else {
        L1 <- sum(diag(Sigma))
    }
    LOD = n.ind/2*log10(exp(1))*(L0 - L1)

    dfM = p * length(terms)
    dfE = p * n.ind - dfM - p
    dfT = p * ( n.ind - 1 )

    SSE = L1
    SST = determinant(crossprod(Y))$modulus
    SSM = SST - SSE


    MSM = SSM / 4
    MSE = SSE / dfE

    result.full <- data.frame( df =c(dfM, dfE,dfT), SS = c(SSM,SSE,SST), MS = c(MSM, MSE,NA), LOD =c( LOD, NA,NA) )
    rownames(result.full) <- c("Model","Error","Total")


    result.drop <- NULL;
    if(length(terms) >= 2) {

        for( drop in terms) {

            X <- NULL
            X <- cbind(rep(1,n.ind))

            for (j in 1:length(terms)) {
                Xadd <- NULL;
                if (length(grep(":", terms[j])) > 0) {
                    if( strsplit(terms[j], " *: *")[[1]][1] %in% names(mark) && strsplit(terms[j], " *: *")[[1]][2] %in% names(mark) && strsplit(terms[j], " *: *")[[1]][1] != drop && strsplit(terms[j], " *: *")[[1]][2] != drop && terms[j] != drop ) {

                        temp <- strsplit(terms[j], " *: *")[[1]]

                        wh1 <- which( strsplit(terms[j], " *: *")[[1]][1] == names(mark))
                        wh2 <- which( strsplit(terms[j], " *: *")[[1]][2] == names(mark))
                        qtt1 <- mark[ wh1 ]
                        qtt2 <- mark[ wh2 ]
                        i1 <- (qtl$chr[ wh1 ])
                        i2 <- (qtl$chr[ wh2 ])

                        Xadd <- cross$geno[[i1]]$prob[ ,qtt1, 1 ] *
                            cross$geno[[i2]]$prob[ ,qtt2, 1 ]

                    }
                } else {
                    if( terms[j] %in% names(mark) && terms[j] != drop ) {
                        wh <- which( terms[j] == names(mark))
                        qtt <- mark[ wh ]
                        i = (qtl$chr[ wh ])

                        Xadd <- cross$geno[[i]]$prob[ ,qtt, 1 ]
                    }
                }

                X <- cbind(X, Xadd)
            }
            E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
            Sigma <- crossprod(E)
            if( method == "hk") {
                L1 <- determinant(Sigma)$modulus
            } else {
                L1 <- sum(diag(Sigma))
            }
            LOD = n.ind/2*log10(exp(1))*(L0 - L1)
            result.drop <- c(result.drop, result.full[1,4] - LOD)
        }
        names(result.drop) <- terms
    }

    return(list( result.full = result.full, result.drop = result.drop, lod = result.full[1,4]) )

}































#refff <- refineqtlM(ril, Y, qtl, formula=formula)


refineqtlM <- function (cross, Y, qtl, chr, pos, qtl.name, formula,
    verbose = TRUE, maxit = 10, incl.markers = TRUE,
    tol = 1e-04, maxit.fitqtl = 1000, method=c("hk","f"), pheno.cols=pheno.cols)
{

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)


    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }


    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    if (!missing(qtl) && (!missing(chr) || !missing(pos) || !missing(qtl.name)))
        warning("qtl argument is provided, and so chr, pos and qtl.name are ignored.")
    if (missing(qtl) && (missing(chr) || missing(pos)))
        stop("Provide either qtl or both chr and pos.")
    if (!missing(qtl)) {
        chr <- qtl$chr
        pos <- qtl$pos
    }  else {
        if (missing(qtl.name)) {
            qtl <- makeqtl(cross, chr = chr, pos = pos, what = "prob")
        }
        else {
            qtl <- makeqtl(cross, chr = chr, pos = pos,
                           qtl.name = qtl.name, what = "prob")
        }
    }

    if (!all(chr %in% names(cross$geno)))
        stop("Chr ", paste(unique(chr[!(chr %in% cross$geno)]),
                           sep = " "), " not found in cross.")

    if (verbose > 1)
        scanqtl.verbose <- TRUE
    else scanqtl.verbose <- FALSE

    cross <- subset(cross, chr = as.character(unique(chr)))

    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what = "prob")
    }
    map <- attr(qtl, "map")

    if (is.null(map))
        stop("Input qtl object should contain the genetic map.")

    mind <- min(sapply(map, function(a) {
        if (is.matrix(a)) a <- a[1, ]
        min(diff(a))
    }))/2
    if (mind <= 0)
        mind <- 1e-06

    if (missing(formula)) {
        formula <- paste("y ~", paste(qtl$altname, collapse = "+"))
        formula <- as.formula(formula)
    }
    formula <- qtl:::checkformula(formula, qtl$altname,NULL)

    tovary <- sort(qtl:::parseformula(formula, qtl$altname, NULL)$idx.qtl)
    if (length(tovary) != qtl$n.qtl)
        reducedqtl <- qtl:::dropfromqtl(qtl, index = (1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl
    if (any(1:length(tovary) != tovary)) {
        tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
        terms <- strsplit(tempform, " *\\+ *")[[1]]
        for (j in seq(along = terms)) {
            if (length(grep(":", terms[j])) > 0) {
                temp <- strsplit(terms[j], " *: *")[[1]]
                for (k in seq(along = temp)) {
                    g <- grep("^[Qq][0-9]+$", temp[k])
                    if (length(g) > 0) {
                        num <- as.numeric(substr(temp[k], 2, nchar(temp[k])))
                        temp[k] <- paste("Q", which(tovary == num),
                                         sep = "")
                  }
                }
                terms[j] <- paste(temp, collapse = ":")
            }
            else {
                g <- grep("^[Qq][0-9]+$", terms[j])
                if (length(g) > 0) {
                    num <- as.numeric(substr(terms[j], 2, nchar(terms[j])))
                    terms[j] <- paste("Q", which(tovary == num),
                                      sep = "")
                }
            }
        }
        formula <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    }
    curpos <- pos[tovary]
    chrnam <- chr[tovary]
    if (verbose)
        cat("pos:", curpos, "\n")
    converged <- FALSE
    oldo <- NULL
    lc <- length(chrnam)
    lastout <- vector("list", length(curpos))
    names(lastout) <- qtl$name[tovary]
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    for (i in 1:maxit) {

        basefit <- getlodM(cross=cross, Y=Y, formula=formula, qtl = qtl, method=method, pheno.cols=pheno.cols)

        if (i == 1) {
            origlod <- curlod <- thisitlod <- as.numeric(basefit)
            origpos <- curpos
        }
        if (verbose)
            cat("Iteration", i, "\n")
        o <- sample(lc)
        if (!is.null(oldo))
            while (o[1] != oldo[lc]) o <- sample(lc)
        oldo <- o
        newpos <- curpos
        for (j in o) {
            otherchr <- chrnam[-j]
            otherpos <- newpos[-j]
            thispos <- as.list(newpos)
            if (any(otherchr == chrnam[j])) {
                linkedpos <- otherpos[otherchr == chr[j]]
                if (any(linkedpos < newpos[j]))
                    low <- max(linkedpos[linkedpos < newpos[j]])
                else low <- -Inf
                if (any(linkedpos > newpos[j]))
                    high <- min(linkedpos[linkedpos > newpos[j]])
                else high <- Inf
                thispos[[j]] <- c(low, high)
            }
            else thispos[[j]] <- c(-Inf, Inf)
            out <- scanqtlM(cross = cross, Y,
                            chr = chrnam, pos = thispos, formula = formula,
                            incl.markers = incl.markers,
                            verbose = scanqtl.verbose, tol = tol, maxit = maxit.fitqtl)
            lastout[[j]] <- out
            newpos[j] <- as.numeric(strsplit(names(out)[out ==
                                                        max(out)], "@")[[1]][2])
            if (verbose) {
                cat(" Q", j, " pos: ", curpos[j], " -> ", newpos[j],
                    "\n", sep = "")
                cat("    LOD increase: ", round(max(out) - curlod,
                                                3), "\n")
            }
            curlod <- max(out)
        }
        if (verbose) {
            cat("all pos:", curpos, "->", newpos, "\n")
            cat("LOD increase at this iteration: ", round(curlod -
                                                          thisitlod, 3), "\n")
        }
        thisitlod <- curlod
        if (max(abs(curpos - newpos)) < mind) {
            converged <- TRUE
            break
        }
        curpos <- newpos
        reducedqtl <- replaceqtl(cross, reducedqtl, seq(length(curpos)),
                                 reducedqtl$chr, curpos, reducedqtl$name)
    }
    if (verbose) {
        cat("overall pos:", origpos, "->", newpos, "\n")
        cat("LOD increase overall: ", round(curlod - origlod,
                                            3), "\n")
    }
    if (!converged)
        warning("Didn't converge.")
    g <- grep("^.+@[0-9\\.]+$", qtl$name)
    if (length(g) == length(qtl$name))
        thenames <- NULL
    else thenames <- qtl$name

    for (j in seq(along = tovary)) qtl <- replaceqtl(cross, qtl,
                  tovary[j], chrnam[j], newpos[j])
    if (!is.null(thenames))
        qtl$name <- thenames

    if ("pLOD" %in% names(attributes(qtl)) && curlod > origlod)
        attr(qtl, "pLOD") <- attr(qtl, "pLOD") + curlod - origlod
    qtl
}






#out <- scanqtlM(cross = cross, Y, chr = chrnam, pos = thispos, formula = formula, incl.markers = incl.markers, verbose = scanqtl.verbose, tol = tol, maxit = maxit.fitqtl)

scanqtlM <- function (cross, Y, chr, pos, formula, incl.markers = FALSE,
                      verbose = TRUE, tol = 1e-04, maxit = 1000, method=c("hk","f"), pheno.cols)
{

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }


    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if (!("prob" %in% names(cross$geno[[1]]))) {
        stop("You need to first run calc.genoprob.")
    }

    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
        attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
        stepwidth.var <- TRUE
        incl.markers <- TRUE
    } else stepwidth.var <- FALSE

    type <- class(cross)[1]
    chrtype <- sapply(cross$geno, class)
    if (length(chr) != length(pos))
        stop("Input chr and pos must have the same length")
    ichr <- match(chr, names(cross$geno))
    if (any(is.na(ichr)))
        stop("There's no chromosome number ", chr[is.na(ichr)],
            " in input cross object")
    n.qtl <- length(chr)
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                sep = "+")
        formula <- as.formula(formula)
    }
    else {
        formula.str <- qtl:::deparseQTLformula(formula)
        for (i in 1:n.qtl) {
            qtl.term <- paste("Q", i, sep = "")
            if (length(grep(qtl.term, formula.str, ignore.case = TRUE)) ==
                0)
                formula.str <- paste(formula.str, qtl.term, sep = "+")
        }
        formula <- as.formula(formula.str)
    }
    formula <- qtl:::checkformula(formula, paste("Q", 1:length(chr),
        sep = ""), NULL)


    sexpgm <- getsex(cross)
    idx.varied <- NULL
    indices <- pos

##

    for (i in 1:length(pos)) {
        l <- length(pos[[i]])
        if (l >= 2) {
            if (l > 2) {
                msg <- "There are more than two elements in "
                msg <- paste(msg, i, "th input pos.")
                msg <- paste(msg, "The first two are taken as starting and ending position.")
                warning(msg)
            }
            idx.varied <- c(idx.varied, i)

            if ("map" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                map <- attr(cross$geno[[ichr[i]]]$prob, "map")
            else {
                stp <- attr(cross$geno[[ichr[i]]]$prob, "step")
                oe <- attr(cross$geno[[ichr[i]]]$prob, "off.end")
                if ("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                    stpw <- attr(cross$geno[[ichr[i]]]$prob,
                                 "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[ichr[i]]]$map,
                                  stp, oe, stpw)

            }
            if (is.matrix(map))
                map <- map[1, ]
            indices[[i]] <- seq(along = map)
            step <- attr(cross$geno[[ichr[i]]]$prob, "step")
            if (!incl.markers && step > 0) {
                eq.sp.pos <- seq(min(map), max(map), by = step)
                wh.eq.pos <- match(eq.sp.pos, map)
                map <- map[wh.eq.pos]
                indices[[i]] <- indices[[i]][wh.eq.pos]
            }
            start <- pos[[i]][1]
            end <- pos[[i]][2]
            tmp <- which((map - start) <= 0)
            if (length(tmp) != 0)
                start <- map[max(tmp)]
            tmp <- which((end - map) <= 0)
            if (length(tmp) != 0)
                end <- map[min(tmp)]
            pos[[i]] <- as.vector(map[(map >= start) & (map <=
                end)])
            indices[[i]] <- indices[[i]][(map >= start) & (map <=
                end)]
        }
    }
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    n.idx.varied <- length(idx.varied)
    n.loop <- 1
    if (n.idx.varied != 0) {
        idx.pos <- rep(0, n.idx.varied)
        l.varied <- NULL
        for (i in 1:n.idx.varied) {
            l.varied[i] <- length(pos[[idx.varied[i]]])
            n.loop <- n.loop * l.varied[i]
        }
        result <- array(rep(0, n.loop), rev(l.varied))
    }  else {
        qtl <- makeqtl(cross, chr = chr, pos = unlist(pos), what = "prob")

        result <- getlodM(cross, Y, formula=formula, qtl = qtl, tol = tol, method = method, pheno.cols=pheno.cols)



        result <- as.numeric(result)
        names(result) <- "LOD"
        class(result) <- "scanqtl"
        attr(result, "method") <- "hk"
        attr(result, "formula") <- qtl:::deparseQTLformula(formula)
        return(result)
    }
    if (verbose) {
        cat(" ", n.loop, "models to fit\n")
        n.prnt <- floor(n.loop/20)
        if (n.prnt < 1)
            n.prnt <- 1
    }
    current.pos <- NULL
    for (i in 1:n.loop) {
        remain <- i
        if (n.idx.varied > 1) {
            for (j in 1:(n.idx.varied - 1)) {
                ns <- 1
                for (k in (j + 1):n.idx.varied) ns <- ns * length(pos[[idx.varied[k]]])
                idx.pos[j] <- floor(remain/ns) + 1
                remain <- remain - (idx.pos[j] - 1) * ns
                if (remain == 0) {
                  idx.pos[j] <- idx.pos[j] - 1
                  remain <- remain + ns
                }
            }
        }
        idx.pos[n.idx.varied] <- remain
        pos.tmp <- NULL
        for (j in 1:length(pos)) {
            if (j %in% idx.varied) {
                idx.tmp <- which(j == idx.varied)
                pos.tmp <- c(pos.tmp, pos[[j]][idx.pos[idx.tmp]])
            }
            else pos.tmp <- c(pos.tmp, pos[[j]])
        }
        if (is.null(current.pos)) {
            qtl.obj <- makeqtl(cross, chr, pos.tmp, what = "prob")
            current.pos <- pos.tmp
        }
        else {
            thew <- rep(NA, length(pos.tmp))
            for (kk in seq(along = pos.tmp)) {
                if (pos.tmp[kk] != current.pos[kk]) {
                  u <- abs(pos.tmp[kk] - pos[[kk]])
                  w <- indices[[kk]][u == min(u)]
                  if (length(w) > 1) {
                    warning("Confused about QTL positions.  You should probably run jittermap to ensure that no two markers conincide.")
                    w <- sample(w, 1)
                }
                  qtl.obj$prob[[kk]] <- cross$geno[[ichr[kk]]]$prob[, w, ]
                  thew[kk] <- w
                  if (chrtype[ichr[kk]] == "X" && (type == "bc" || type == "f2")) {

                      temp <- qtl.obj$prob[[kk]]
                      temp <- array(temp, dim = c(nrow(temp), 1, ncol(temp)))
                      dimnames(temp) <- list(NULL, "loc", 1:ncol(qtl.obj$prob[[kk]]))
                      qtl.obj$prob[[kk]] <- qtl:::reviseXdata(type,
                      "full", sexpgm, prob = temp, cross.attr = attributes(cross))[, 1, ]
                  }
                  current.pos[kk] <- pos.tmp[kk]
              }
            }
        }

        qtl.obj <- makeqtl(cross, chr, current.pos, what = "prob")
        fit <- getlodM(cross, Y, formula=formula, qtl = qtl.obj, tol=tol, method=method, pheno.cols=pheno.cols)

        if (verbose && ((i - 1)%%n.prnt) == 0)
            cat("    ", i, "/", n.loop, "\n")
        result[i] <- as.numeric(fit)
    }
    dnames <- list(NULL)
    for (i in 1:n.idx.varied) {
        i.chr <- chr[idx.varied[n.idx.varied - i + 1]]
        i.pos <- pos[[idx.varied[n.idx.varied - i + 1]]]
        dnames[[i]] <- paste(paste("Chr", i.chr, sep = ""), i.pos,
                             sep = "@")
    }
    dimnames(result) <- dnames
    class(result) <- "scanqtl"
    attr(result, "method") <- "hk"
    attr(result, "formula") <- qtl:::deparseQTLformula(formula)
    result
}




#out <- addqtlM(cross, Y = Y, qtl = qtl,
#                  formula = thisformula,
#                  incl.markers = incl.markers, verbose = verbose.scan)

#
#out <- addqtlM(cross, Y = Y, qtl = qtl,  formula = firstformula)

#out2 <- addqtl(cross, pheno.col = 1, qtl = qtl,
#              covar = NULL, formula = thisformula, method = method,
#              incl.markers = incl.markers, verbose = verbose.scan)

#out <- addqtl(cross, pheno.col = 1, qtl = qtl,
#              covar = NULL, formula = firstformula, method = method,
#              incl.markers = incl.markers, verbose = verbose.scan)




addqtlM <- function (cross, Y, chr, qtl, formula,
    incl.markers = TRUE, verbose = FALSE, tol = 1e-04, maxit = 1000, method = c("hk","f"), pheno.cols)
{


    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")
    if (!("qtl" %in% class(qtl)))
        stop("The qtl argument must be an object of class \"qtl\".")
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if (verbose > 1) {
        verbose <- TRUE
        verbose.scanqtl <- TRUE
    } else verbose.scanqtl <- FALSE

    n.qtl <- qtl$n.qtl
    qtlchr <- qtl$chr
    qtlpos <- qtl$pos

    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }
    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
        attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
        stepwidth.var <- TRUE
        incl.markers <- TRUE
    } else stepwidth.var <- FALSE


    if (missing(chr)) {
        chr <- names(cross$geno)
    } else chr <- qtl:::matchchr(chr, names(cross$geno))


    n.covar <- 0
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        newformula <- as.formula(paste(formula, "+Q", n.qtl +
            1, sep = ""))
        formula <- as.formula(formula)
    }
    else {
        newqtl <- paste("Q", n.qtl + 1, sep = "")
        formula <- qtl:::checkformula(formula, c(qtl$altname, newqtl),
                                NULL)
        theterms <- rownames(attr(terms(formula), "factors"))
        g <- grep(paste("^[Qq]", n.qtl + 1, "$", sep = ""), theterms)
        if (length(g) == 0) {
            newformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                                           "+ Q", n.qtl + 1, sep = ""))
        }     else {
            newformula <- formula
            theterms <- colnames(attr(terms(formula), "factors"))
            g <- unique(c(grep(paste("^[Qq]", n.qtl + 1, "$",
                sep = ""), theterms), grep(paste("^[Qq]", n.qtl +
                1, " *:", sep = ""), theterms), grep(paste(": +[Qq]",
                n.qtl + 1, " *:", sep = ""), theterms), grep(paste(": +[Qq]",
                n.qtl + 1, "$", sep = ""), theterms)))
            if (length(g) > 0) {
                theterms <- theterms[-g]
                formula <- as.formula(paste("y ~ ", paste(theterms,
                  collapse = " + "), sep = ""))
            }
        }
    }
    thefactors <- rownames(attr(terms(formula), "factors"))
    todrop <- NULL
    for (i in 1:n.qtl) {
        if (length(grep(paste("^[Qq]", i, "$", sep = ""), thefactors)) ==
            0)
            todrop <- c(todrop, i)
    }
    if (length(todrop) > 0) {
        newqtlnum <- n.qtl + 1
        notdropped <- (1:n.qtl)[-todrop]
        newnum <- 1:length(notdropped)
        qtl <- dropfromqtl(qtl, index = todrop)
        qtlchr <- qtlchr[-todrop]
        qtlpos <- qtlpos[-todrop]
        n.qtl <- n.qtl - length(todrop)
        revnewqtlnum <- n.qtl + 1
        formula <- qtl:::reviseqtlnuminformula(formula, notdropped,
            newnum)
        newformula <- qtl:::reviseqtlnuminformula(newformula, c(notdropped,
            newqtlnum), c(newnum, revnewqtlnum))
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    lod0 <- as.numeric(getlodM(cross, Y, qtl = qtl, formula = formula, tol= tol, method=method, pheno.cols=pheno.cols))



#lod0 <- qtl:::fitqtlengine(pheno = pheno[,1], qtl = qtl, covar = NULL,
#        formula = formula, method = "hk", model = "normal", dropone = FALSE,
#        get.ests = FALSE, run.checks = FALSE, cross.attr = cross.attr,
#        sexpgm = sexpgm, tol = tol, maxit = maxit)$result.full[1, 4]



    results <- NULL
    for (i in chr) {
        if (verbose)
            cat("Scanning chr", i, "\n")
        thechr <- c(qtlchr, i)
        thepos <- c(as.list(qtlpos), list(c(-Inf, Inf)))
        sqout <- scanqtlM(cross, Y, chr = thechr,
            pos = thepos, formula = newformula,
            incl.markers = incl.markers,
            verbose = verbose.scanqtl, tol = tol, maxit = maxit, method=method, pheno.cols=pheno.cols)

        if ("map" %in% names(attributes(cross$geno[[i]]$prob)))
            map <- attr(cross$geno[[i]]$prob, "map")
        else {
            stp <- attr(cross$geno[[i]]$prob, "step")
            oe <- attr(cross$geno[[i]]$prob, "off.end")
            if ("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
            else stpw <- "fixed"
            map <- create.map(cross$geno[[i]]$map, stp, oe,
                              stpw)
        }

        if (is.matrix(map))
            map <- map[1, ]
        step <- attr(cross$geno[[i]]$prob, "step")
        if (!incl.markers && step > 0) {
            eq.sp.pos <- seq(min(map), max(map), by = step)
            wh.eq.pos <- match(eq.sp.pos, map)
            map <- map[wh.eq.pos]
        }
        w <- names(map)
        o <- grep("^loc-*[0-9]+", w)
        if (length(o) > 0)
            w[o] <- paste("c", i, ".", w[o], sep = "")
        z <- data.frame(lod = as.numeric(sqout) - lod0, stringsAsFactors = TRUE)
        z <- cbind(chr = rep(i, length(map)), pos = as.numeric(map),
            z)
        rownames(z) <- w
        results <- rbind(results, z)
    }
    class(results) <- c("scanone", "data.frame")
    attr(results, "method") <- "hk"
    attr(results, "formula") <- qtl:::deparseQTLformula(newformula)
    results
}







addintM <- function (cross, Y, qtl, formula, qtl.only = FALSE, verbose = TRUE,
    pvalues = TRUE, simple = FALSE, tol = 1e-04, maxit = 1000, method=method, pheno.cols=pheno.cols)
{

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")
    if (!("qtl" %in% class(qtl)))
        stop("The qtl argument must be an object of class \"qtl\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)
    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }

    n.qtl <- qtl$n.qtl
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                sep = "+")
        if (n.covar) {
            tmp.C <- colnames(covar)
            for (i in 1:n.covar) formula <- paste(formula, tmp.C[i],
                sep = "+")
        }
        formula <- as.formula(formula)
    }
    formula <- qtl:::checkformula(formula, qtl$altname, NULL)
    factors <- attr(terms(formula), "factors")
    if (sum(factors[1, ]) == 0)
        factors <- factors[-1, ]
    fn <- fn.alt <- rownames(factors)
    qan <- qtl$altname
    qn <- qtl$name
    m <- match(fn, qan)
    fn.alt[!is.na(m)] <- qn[m[!is.na(m)]]
    int2test <- int2test.alt <- NULL
    for (i in 1:(nrow(factors) - 1)) {
        for (j in (i + 1):nrow(factors)) {
            temp <- rep(0, nrow(factors))
            temp[c(i, j)] <- 1
            if (!any(apply(factors, 2, function(a, b) all(a ==
                b), temp))) {
                int2test <- c(int2test, paste(fn[i], fn[j], sep = ":"))
                int2test.alt <- c(int2test.alt, paste(fn.alt[i],
                  fn.alt[j], sep = ":"))
            }
        }
    }
    if (qtl.only && length(int2test) > 0) {
        z <- matrix(unlist(strsplit(int2test, ":")), ncol = 2,
            byrow = TRUE)
        wh <- apply(z, 1, function(a) length(grep("^[Qq][0-9]+$",
            a)))
        int2test <- int2test[wh == 2]
        int2test.alt <- int2test.alt[wh == 2]
    }
    n2test <- length(int2test)
    if (n2test == 0) {
        if (verbose)
            cat("No pairwise interactions to add.\n")
        return(NULL)
    }
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

###
    thefit0 <- fitqtlengineM(cross = cross, Y = Y, qtl = qtl,
        formula = formula, tol = tol, method=method, pheno.cols=pheno.cols)

    results <- matrix(ncol = 7, nrow = n2test)
    dimnames(results) <- list(int2test.alt, c("df", "Type III SS",
        "LOD", "%var", "F value", "Pvalue(Chi2)", "Pvalue(F)"))
    for (k in seq(along = int2test)) {


        thefit1 <- fitqtlengineM(cross = cross, Y = Y, qtl = qtl,
            formula = as.formula(paste(qtl:::deparseQTLformula(formula),
                int2test[k], sep = "+")), tol = tol, method=method, pheno.cols=pheno.cols)

        results[k, 1] <- thefit1$result.full[1, 1] - thefit0$result.full[1,
            1]
        results[k, 2] <- thefit1$result.full[1, 2] - thefit0$result.full[1,
            2]
        results[k, 3] <- thefit1$result.full[1, 4] - thefit0$result.full[1,
            4]
        results[k, 4] <- 100 * (1 - 10^(-2 * thefit1$result.full[1,
            4]/qtl$n.ind)) - 100 * (1 - 10^(-2 * thefit0$result.full[1,
            4]/qtl$n.ind))
        results[k, 5] <- (results[k, 2]/results[k, 1])/thefit1$result.full[2,
            3]
        results[k, 6] <- pchisq(results[k, 3] * 2 * log(10),
            results[k, 1], lower.tail = FALSE)
        results[k, 7] <- pf(results[k, 5], results[k, 1], thefit1$result.full[3,
            1], lower.tail = FALSE)
    }
    results <- as.data.frame(results, stringsAsFactors = TRUE)
    class(results) <- c("addint", "data.frame")
    attr(results, "method") <- "hk"
    attr(results, "model") <- "normal"
    attr(results, "formula") <- qtl:::deparseQTLformula(formula)
    if (simple)
        pvalues <- FALSE
    attr(results, "pvalues") <- pvalues
    attr(results, "simple") <- simple
    results
}








#outt <- stepwiseqtlM(ril, Y=Y, max.qtl = 7, penalties=c(3,2,2) )



stepwiseqtlM <- function (cross, chr, Y, qtl, formula, max.qtl = 10,
    incl.markers = TRUE, refine.locations = TRUE,
    penalties,  additive.only = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000, method=method, pheno.cols=pheno.cols)
{


    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }
    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")
    if (!missing(chr))
        cross <- subset(cross, chr)

    if (!missing(qtl)) {
        if (!("qtl" %in% class(qtl)))
            stop("The qtl argument must be an object of class \"qtl\".")
        m <- is.na(match(qtl$chr, names(cross$geno)))
        if (any(m)) {
            wh <- qtl$chr[m]
            if (length(wh) > 1)
                stop("Chromosomes ", paste(wh, collapse = ", "),
                  " (in QTL object) not in cross object.")
            else stop("Chromosome ", wh, " (in QTL object) not in cross object.")
        }
        if (missing(formula)) {
            if (!is.null(covar))
                formula <- paste("y ~ ", paste(names(covar),
                  collapse = "+"), "+")
            else formula <- "y ~ "
            formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr),
                sep = ""), collapse = "+"))
        }
        else {
            temp <- qtl:::checkStepwiseqtlStart(qtl, formula, covar)
            qtl <- temp$qtl
            formula <- temp$formula
        }
        startatnull <- FALSE
    }
    else {
        if (!missing(formula))
            warning("formula ignored if qtl is not provided.")
        startatnull <- TRUE
    }
    if (!startatnull)
        qtl$name <- qtl$altname

    qtlmethod <- "prob"
    if (!missing(qtl) && qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }

    if (!startatnull) {
        if (method == "hk" && !("prob" %in% names(qtl)))
            stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }

    if (max.qtl < 1)
        stop("Need max.qtl > 0 if we are to scan for qtl")
    lod0 <- 0
    if (startatnull)
        firstformula <- y ~ Q1
    else firstformula <- formula

    cross.type <- class(cross)[1]
    if (missing(penalties)) {
        if (cross.type == "f2") {
            penalties <- c(3.52, 4.28, 2.69)
        }
        else if (cross.type == "bc") {
            penalties <- c(2.69, 2.62, 1.19)
        }
        else stop("No default penalties available for cross type ",
            cross.type)
    }
    else if (length(penalties) != 3) {
        if (length(penalties) == 1) {
            stop("You must include a penalty for interaction terms.")
        }
        else {
            if (length(penalties) == 2)
                penalties <- penalties[c(1, 2, 2)]
            else {
                warning("penalties should have length 3")
                penalties <- penalties[1:3]
            }
        }
    }
    if (verbose > 2)
        verbose.scan <- TRUE
    else verbose.scan <- FALSE

    curbest <- NULL
    curbestplod <- 0
    if (verbose)
        cat(" -Initial scan\n")
    if (startatnull) {
        {
            out <- scanoneM(cross, Y)
            lod <- max(out[, 3], na.rm = TRUE)
            curplod <- calc.plod(lod, c(1, 0, 0), penalties = penalties)
            wh <- which(!is.na(out[, 3]) & out[, 3] == lod)
            if (length(wh) > 1)
                wh <- sample(wh, 1)
            qtl <- makeqtl(cross, as.character(out[wh, 1]), out[wh,
                2], "Q1", what = qtlmethod)
            formula <- firstformula
            n.qtl <- 1
        }
        if (max.qtl !=1 ) {

            out <- addqtlM(cross, Y = Y, qtl = qtl,
                          formula = firstformula,
            incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)


            curlod <- max(out[, 3], na.rm = TRUE)
            wh <- which(!is.na(out[, 3]) & out[, 3] == curlod)
            if (length(wh) > 1)
                wh <- sample(wh, 1)
            curqtl <- addtoqtl(cross, qtl, as.character(out[wh, 1]),
                               out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
            curformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                 "+Q", n.qtl + 1, sep = ""))
            curlod <- curlod + lod
            curplod <- calc.plod(curlod, qtl:::countqtlterms(curformula,
                 ignore.covar = TRUE), penalties = penalties)

#            if (verbose)
#                cat("        plod =", curplod, "\n")
            curnqtl <- n.qtl + 1

            if (!additive.only) {

#
                thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                               "+Q", n.qtl + 1, "+Q", 1, ":Q", n.qtl + 1,
                               sep = ""))
                out <- addqtlM(cross, Y = Y, qtl = qtl,
                  formula = thisformula,
                  incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)

#


                thislod <- max(out[, 3], na.rm = TRUE)
                wh <- which(!is.na(out[, 3]) & out[, 3] == thislod)
                if (length(wh) > 1)
                  wh <- sample(wh, 1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,
                  1]), out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
                thislod <- thislod + lod
                thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                  ignore.covar = TRUE), penalties = penalties)
#                if (verbose)
#                  cat("        plod =", thisplod, "\n")
                if (thisplod > curplod) {
                  curformula <- thisformula
                  curplod <- thisplod
                  curlod <- thislod
                  curqtl <- thisqtl
                  curnqtl <- n.qtl + 1
                }
            }



            n.qtl <- curnqtl
            qtl <- curqtl
            formula <- curformula
            lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
                                      tol = tol, method=method, pheno.cols=pheno.cols)) - lod0

        }
    }  else {
        if (verbose)
            cat(" ---Starting at a model with", length(qtl$chr),
                "QTL\n")
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
            }
            qtl <- rqtl
        }
        lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
            tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        attr(qtl, "pLOD") <- curplod
        n.qtl <- length(qtl$chr)
    }
    attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod
    if (curplod > 0) {
        curbest <- qtl
        curbestplod <- curplod
        if (verbose)
            cat("** new best ** (pLOD increased by ", round(curplod,
                4), ")\n", sep = "")
    }

    if (verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            qtl:::deparseQTLformula(formula), "\n")
    if (verbose > 1)
        cat("         qtl:", paste(qtl$chr, round(qtl$pos, 1),
            sep = "@"), "\n")
    i <- 0
    while (n.qtl < max.qtl) {
        i <- i + 1
        if (verbose) {
            cat(" -Step", i, "\n")
            cat(" ---Scanning for additive qtl\n")
        }
        out <- addqtlM(cross, Y = Y, qtl = qtl,
            formula = formula, incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
        curlod <- max(out[, 3], na.rm = TRUE)
        wh <- which(!is.na(out[, 3]) & out[, 3] == curlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        curqtl <- addtoqtl(cross, qtl, as.character(out[wh, 1]),
            out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
        curformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
            "+Q", n.qtl + 1, sep = ""))
        curlod <- curlod + lod
        curplod <- calc.plod(curlod, qtl:::countqtlterms(curformula,
            ignore.covar = TRUE), penalties = penalties)
        if (verbose)
            cat("        plod =", curplod, "\n")
        curnqtl <- n.qtl + 1
        if (!additive.only) {
            for (j in 1:n.qtl) {
                if (verbose)
                  cat(" ---Scanning for QTL interacting with Q",
                    j, "\n", sep = "")
                thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                  "+Q", n.qtl + 1, "+Q", j, ":Q", n.qtl + 1,
                  sep = ""))
                out <- addqtlM(cross, Y = Y, qtl = qtl,
                  formula = thisformula,
                  incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
                thislod <- max(out[, 3], na.rm = TRUE)
                wh <- which(!is.na(out[, 3]) & out[, 3] == thislod)
                if (length(wh) > 1)
                  wh <- sample(wh, 1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,
                  1]), out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
                thislod <- thislod + lod
                thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                  ignore.covar = TRUE), penalties = penalties)
                if (verbose)
                  cat("        plod =", thisplod, "\n")
                if (thisplod > curplod) {
                  curformula <- thisformula
                  curplod <- thisplod
                  curlod <- thislod
                  curqtl <- thisqtl
                  curnqtl <- n.qtl + 1
                }
            }
            if (n.qtl > 1) {
                if (verbose)
                  cat(" ---Look for additional interactions\n")



                temp <- addintM(cross, Y = Y, qtl,
                  formula = formula, qtl.only = TRUE,
                  verbose = verbose.scan, method=method, pheno.cols=pheno.cols)






                if (!is.null(temp)) {
                  thislod <- max(temp[, 3], na.rm = TRUE)
                  wh <- which(!is.na(temp[, 3]) & temp[, 3] ==
                    thislod)
                  if (length(wh) > 1)
                    wh <- sample(wh, 1)
                  thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                    "+", rownames(temp)[wh]))
                  thislod <- thislod + lod
                  thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                    ignore.covar = TRUE), penalties = penalties)
                  if (verbose)
                    cat("        plod =", thisplod, "\n")
                  if (thisplod > curplod) {
                    curformula <- thisformula
                    curplod <- thisplod
                    curlod <- thislod
                    curqtl <- qtl
                    curnqtl <- n.qtl
                  }
                }
            }

        }
        qtl <- curqtl
        n.qtl <- curnqtl
        attr(qtl, "formula") <- qtl:::deparseQTLformula(curformula)
        attr(qtl, "pLOD") <- curplod
        formula <- curformula
        lod <- curlod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
                qtl <- rqtl
                lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
                  tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
                curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                  ignore.covar = TRUE), penalties = penalties)
                attr(qtl, "pLOD") <- curplod
            }
        }
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,
                1), sep = "@"), "\n")
        if (curplod > curbestplod) {
            if (verbose)
                cat("** new best ** (pLOD increased by ", round(curplod -
                  curbestplod, 4), ")\n", sep = "")
            curbest <- qtl
            curbestplod <- curplod
        }
        if (n.qtl >= max.qtl)
            break
    }
    if (verbose)
        cat(" -Starting backward deletion\n")
    while (n.qtl > 1) {
        i <- i + 1

###
        out <- fitqtlM(cross, Y = Y, qtl = qtl, formula = formula, method=method, pheno.cols=pheno.cols)$result.drop


        rn <- names(out)
        wh <- c(grep("^[Qq][0-9]+$", rn), grep("^[Qq][0-9]+:[Qq][0-9]+$",
            rn))
        thelod <- out[wh]
        minlod <- min(thelod, na.rm = TRUE)
        wh <- which(!is.na(thelod) & thelod == minlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        lod <- lod - minlod
        todrop <- names(out)[wh]
        if (verbose)
            cat(" ---Dropping", todrop, "\n")
        if (length(grep(":", todrop)) > 0) {
            theterms <- attr(terms(formula), "factors")
            wh <- colnames(theterms) == todrop
            if (!any(wh))
                stop("Confusion about what interation to drop!")
            theterms <- colnames(theterms)[!wh]
            formula <- as.formula(paste("y~", paste(theterms,
                collapse = "+")))
        }
        else {
            numtodrop <- as.numeric(substr(todrop, 2, nchar(todrop)))
            theterms <- attr(terms(formula), "factors")
            cn <- colnames(theterms)
            g <- c(grep(paste("^[Qq]", numtodrop, "$", sep = ""),
                cn), grep(paste("^[Qq]", numtodrop, ":", sep = ""),
                cn), grep(paste(":[Qq]", numtodrop, "$", sep = ""),
                cn))
            cn <- cn[-g]
            formula <- as.formula(paste("y~", paste(cn, collapse = "+")))
            if (n.qtl > numtodrop) {
                for (j in (numtodrop + 1):n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                  j, j - 1)
            }
            qtl <- dropfromqtl(qtl, index = numtodrop)
            qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl,
                sep = "")
            n.qtl <- n.qtl - 1
        }
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,
                1), sep = ":"), "\n")
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- curplod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            if (!is.null(qtl)) {


                rqtl <- refineqtlM(cross, Y = Y,
                  qtl = qtl, formula = formula,
                  verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)

                if (any(rqtl$pos != qtl$pos)) {
                  if (verbose)
                    cat(" ---  Moved a bit\n")
                  qtl <- rqtl
                  lod <- as.numeric( getlodM(cross, Y = Y, qtl = qtl,
                    formula = formula, tol = tol, method=method, pheno.cols=pheno.cols) ) - lod0
                  curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                    ignore.covar = TRUE), penalties = penalties)
                  attr(qtl, "pLOD") <- curplod
                }
            }
        }
        if (curplod > curbestplod) {
            if (verbose)
                cat("** new best ** (pLOD increased by ", round(curplod -
                  curbestplod, 4), ")\n", sep = "")
            curbestplod <- curplod
            curbest <- qtl
        }
    }
    if (!is.null(curbest)) {
        chr <- curbest$chr
        pos <- curbest$pos
        o <- order(factor(chr, levels = names(cross$geno)), pos)
        qtl <- makeqtl(cross, chr[o], pos[o], what = qtlmethod)
        formula <- as.formula(attr(curbest, "formula"))
        if (length(chr) > 1) {
            n.qtl <- length(chr)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                i, n.qtl + i)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                n.qtl + o[i], i)
        }
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- attr(curbest, "pLOD")
        curbest <- qtl
    }
    else {
        curbest <- numeric(0)
        class(curbest) <- "qtl"
        attr(curbest, "pLOD") <- 0
    }
    attr(curbest, "formula") <- qtl:::deparseQTLformula(attr(curbest,
        "formula"), TRUE)
    curbest
}














stepwiseqtlM <- function (cross, chr, Y, qtl, formula, max.qtl = 10,
    incl.markers = TRUE, refine.locations = TRUE,
    penalties,  additive.only = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000, method=method, pheno.cols=pheno.cols)
{
    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")
    if (!missing(chr))
        cross <- subset(cross, chr)

    if (!missing(qtl)) {
        if (!("qtl" %in% class(qtl)))
            stop("The qtl argument must be an object of class \"qtl\".")
        m <- is.na(match(qtl$chr, names(cross$geno)))
        if (any(m)) {
            wh <- qtl$chr[m]
            if (length(wh) > 1)
                stop("Chromosomes ", paste(wh, collapse = ", "),
                  " (in QTL object) not in cross object.")
            else stop("Chromosome ", wh, " (in QTL object) not in cross object.")
        }
        if (missing(formula)) {
            if (!is.null(covar))
                formula <- paste("y ~ ", paste(names(covar),
                  collapse = "+"), "+")
            else formula <- "y ~ "
            formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr),
                sep = ""), collapse = "+"))
        }
        else {
            temp <- qtl:::checkStepwiseqtlStart(qtl, formula, covar)
            qtl <- temp$qtl
            formula <- temp$formula
        }
        startatnull <- FALSE
    }
    else {
        if (!missing(formula))
            warning("formula ignored if qtl is not provided.")
        startatnull <- TRUE
    }
    if (!startatnull)
        qtl$name <- qtl$altname

    qtlmethod <- "prob"
    if (!missing(qtl) && qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }

    if (!startatnull) {
        if (method == "hk" && !("prob" %in% names(qtl)))
            stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }

    if (max.qtl < 1)
        stop("Need max.qtl > 0 if we are to scan for qtl")
    lod0 <- 0
    if (startatnull)
        firstformula <- y ~ Q1
    else firstformula <- formula

    cross.type <- class(cross)[1]
    if (missing(penalties)) {
        if (cross.type == "f2") {
            penalties <- c(3.52, 4.28, 2.69)
        }
        else if (cross.type == "bc") {
            penalties <- c(2.69, 2.62, 1.19)
        }
        else stop("No default penalties available for cross type ",
            cross.type)
    }
    else if (length(penalties) != 3) {
        if (length(penalties) == 1) {
            stop("You must include a penalty for interaction terms.")
        }
        else {
            if (length(penalties) == 2)
                penalties <- penalties[c(1, 2, 2)]
            else {
                warning("penalties should have length 3")
                penalties <- penalties[1:3]
            }
        }
    }
    if (verbose > 2)
        verbose.scan <- TRUE
    else verbose.scan <- FALSE

    curbest <- NULL
    curbestplod <- 0
    if (verbose)
        cat(" -Initial scan\n")
    if (startatnull) {
        {
            out <- scanoneM(cross, Y)
            lod <- max(out[, 3], na.rm = TRUE)
            curplod <- calc.plod(lod, c(1, 0, 0), penalties = penalties)
            wh <- which(!is.na(out[, 3]) & out[, 3] == lod)
            if (length(wh) > 1)
                wh <- sample(wh, 1)
            qtl <- makeqtl(cross, as.character(out[wh, 1]), out[wh,
                2], "Q1", what = qtlmethod)
            formula <- firstformula
            n.qtl <- 1
        }

    }  else {
        if (verbose)
            cat(" ---Starting at a model with", length(qtl$chr),
                "QTL\n")
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
            }
            qtl <- rqtl
        }
        lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
            tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        attr(qtl, "pLOD") <- curplod
        n.qtl <- length(qtl$chr)
    }
    attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod
    if (curplod > 0) {
        curbest <- qtl
        curbestplod <- curplod
        if (verbose)
            cat("** new best ** (pLOD increased by ", round(curplod,
                4), ")\n", sep = "")
    }

    if (verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            qtl:::deparseQTLformula(formula), "\n")
    if (verbose > 1)
        cat("         qtl:", paste(qtl$chr, round(qtl$pos, 1),
            sep = "@"), "\n")
    i <- 0
    while (n.qtl < max.qtl) {
        i <- i + 1
        if (verbose) {
            cat(" -Step", i, "\n")
            cat(" ---Scanning for additive qtl\n")
        }
        out <- addqtlM(cross, Y = Y, qtl = qtl,
            formula = formula, incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
        curlod <- max(out[, 3], na.rm = TRUE)
        wh <- which(!is.na(out[, 3]) & out[, 3] == curlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        curqtl <- addtoqtl(cross, qtl, as.character(out[wh, 1]),
            out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
        curformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
            "+Q", n.qtl + 1, sep = ""))
        curlod <- curlod + lod
        curplod <- calc.plod(curlod, qtl:::countqtlterms(curformula,
            ignore.covar = TRUE), penalties = penalties)
        if (verbose)
            cat("        plod =", curplod, "\n")
        curnqtl <- n.qtl + 1
        if (!additive.only) {
            for (j in 1:n.qtl) {
                if (verbose)
                  cat(" ---Scanning for QTL interacting with Q",
                    j, "\n", sep = "")
                thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                  "+Q", n.qtl + 1, "+Q", j, ":Q", n.qtl + 1,
                  sep = ""))
                out <- addqtlM(cross, Y = Y, qtl = qtl,
                  formula = thisformula,
                  incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
                thislod <- max(out[, 3], na.rm = TRUE)
                wh <- which(!is.na(out[, 3]) & out[, 3] == thislod)
                if (length(wh) > 1)
                  wh <- sample(wh, 1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,
                  1]), out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
                thislod <- thislod + lod
                thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                  ignore.covar = TRUE), penalties = penalties)
                if (verbose)
                  cat("        plod =", thisplod, "\n")
                if (thisplod > curplod) {
                  curformula <- thisformula
                  curplod <- thisplod
                  curlod <- thislod
                  curqtl <- thisqtl
                  curnqtl <- n.qtl + 1
                }
            }
            if (n.qtl > 1) {
                if (verbose)
                  cat(" ---Look for additional interactions\n")



                temp <- addintM(cross, Y = Y, qtl,
                  formula = formula, qtl.only = TRUE,
                  verbose = verbose.scan, method=method, pheno.cols=pheno.cols)






                if (!is.null(temp)) {
                  thislod <- max(temp[, 3], na.rm = TRUE)
                  wh <- which(!is.na(temp[, 3]) & temp[, 3] ==
                    thislod)
                  if (length(wh) > 1)
                    wh <- sample(wh, 1)
                  thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                    "+", rownames(temp)[wh]))
                  thislod <- thislod + lod
                  thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                    ignore.covar = TRUE), penalties = penalties)
                  if (verbose)
                    cat("        plod =", thisplod, "\n")
                  if (thisplod > curplod) {
                    curformula <- thisformula
                    curplod <- thisplod
                    curlod <- thislod
                    curqtl <- qtl
                    curnqtl <- n.qtl
                  }
                }
            }

        }
        qtl <- curqtl
        n.qtl <- curnqtl
        attr(qtl, "formula") <- qtl:::deparseQTLformula(curformula)
        attr(qtl, "pLOD") <- curplod
        formula <- curformula
        lod <- curlod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
                qtl <- rqtl
                lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
                  tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
                curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                  ignore.covar = TRUE), penalties = penalties)
                attr(qtl, "pLOD") <- curplod
            }
        }
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,
                1), sep = "@"), "\n")
        if (curplod > curbestplod) {
            if (verbose)
                cat("** new best ** (pLOD increased by ", round(curplod -
                  curbestplod, 4), ")\n", sep = "")
            curbest <- qtl
            curbestplod <- curplod
        }
        if (n.qtl >= max.qtl)
            break
    }
    if (verbose)
        cat(" -Starting backward deletion\n")
    while (n.qtl > 1) {
        i <- i + 1

###
        out <- fitqtlM(cross, Y = Y, qtl = qtl, formula = formula, method=method, pheno.cols=pheno.cols)$result.drop


        rn <- names(out)
        wh <- c(grep("^[Qq][0-9]+$", rn), grep("^[Qq][0-9]+:[Qq][0-9]+$",
            rn))
        thelod <- out[wh]
        minlod <- min(thelod, na.rm = TRUE)
        wh <- which(!is.na(thelod) & thelod == minlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        lod <- lod - minlod
        todrop <- names(out)[wh]
        if (verbose)
            cat(" ---Dropping", todrop, "\n")
        if (length(grep(":", todrop)) > 0) {
            theterms <- attr(terms(formula), "factors")
            wh <- colnames(theterms) == todrop
            if (!any(wh))
                stop("Confusion about what interation to drop!")
            theterms <- colnames(theterms)[!wh]
            formula <- as.formula(paste("y~", paste(theterms,
                collapse = "+")))
        }
        else {
            numtodrop <- as.numeric(substr(todrop, 2, nchar(todrop)))
            theterms <- attr(terms(formula), "factors")
            cn <- colnames(theterms)
            g <- c(grep(paste("^[Qq]", numtodrop, "$", sep = ""),
                cn), grep(paste("^[Qq]", numtodrop, ":", sep = ""),
                cn), grep(paste(":[Qq]", numtodrop, "$", sep = ""),
                cn))
            cn <- cn[-g]
            formula <- as.formula(paste("y~", paste(cn, collapse = "+")))
            if (n.qtl > numtodrop) {
                for (j in (numtodrop + 1):n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                  j, j - 1)
            }
            qtl <- dropfromqtl(qtl, index = numtodrop)
            qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl,
                sep = "")
            n.qtl <- n.qtl - 1
        }
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,
                1), sep = ":"), "\n")
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- curplod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            if (!is.null(qtl)) {


                rqtl <- refineqtlM(cross, Y = Y,
                  qtl = qtl, formula = formula,
                  verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)

                if (any(rqtl$pos != qtl$pos)) {
                  if (verbose)
                    cat(" ---  Moved a bit\n")
                  qtl <- rqtl
                  lod <- as.numeric( getlodM(cross, Y = Y, qtl = qtl,
                    formula = formula, tol = tol, method=method, pheno.cols=pheno.cols) ) - lod0
                  curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                    ignore.covar = TRUE), penalties = penalties)
                  attr(qtl, "pLOD") <- curplod
                }
            }
        }
        if (curplod > curbestplod) {
            if (verbose)
                cat("** new best ** (pLOD increased by ", round(curplod -
                  curbestplod, 4), ")\n", sep = "")
            curbestplod <- curplod
            curbest <- qtl
        }
    }
    if (!is.null(curbest)) {
        chr <- curbest$chr
        pos <- curbest$pos
        o <- order(factor(chr, levels = names(cross$geno)), pos)
        qtl <- makeqtl(cross, chr[o], pos[o], what = qtlmethod)
        formula <- as.formula(attr(curbest, "formula"))
        if (length(chr) > 1) {
            n.qtl <- length(chr)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                i, n.qtl + i)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                n.qtl + o[i], i)
        }
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- attr(curbest, "pLOD")
        curbest <- qtl
    }
    else {
        curbest <- numeric(0)
        class(curbest) <- "qtl"
        attr(curbest, "pLOD") <- 0
    }
    attr(curbest, "formula") <- qtl:::deparseQTLformula(attr(curbest,
        "formula"), TRUE)
    curbest
}


