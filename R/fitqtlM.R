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
        Y <- as.matrix(cross$pheno[,pheno.cols])
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    E <- matrix(NA, n.ind, p)
    X <- cbind(rep(1,n.ind))
    E <- stats::lm.fit(X, Y, tol=tol)$residuals
    Sigma <- crossprod(E)
    if( method == "hk") {
        L0 <- determinant(Sigma)$modulus
    } else {
        L0 <- sum(diag(Sigma))
    }


    tempform <- strsplit(qtl::deparseQTLformula(formula), " *~ *")[[1]][2]
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
    E <- lm.fit(X, Y, tol=tol)$residuals
    Sigma <- crossprod(E)
    if( method == "hk") {
        L1 <- determinant(Sigma)$modulus
        LOD <- n.ind/2*(L0 - L1)/log(10)
    } else {
        L1 <- sum(diag(Sigma))
        LOD <- n.ind/2*log(L0/L1, 10)
    }



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
            E <- lm.fit(X, Y, tol=tol)$residuals
            Sigma <- crossprod(E)
            if( method == "hk") {
                L1 <- determinant(Sigma)$modulus
                LOD <- n.ind/2*(L0 - L1)/log(10)
            } else {
                L1 <- sum(diag(Sigma))
                LOD <- n.ind/2*log(L0/L1, 10)
            }
            result.drop <- c(result.drop, result.full[1,4] - LOD)
        }
        names(result.drop) <- qtl$name
    }

    return(list( result.full = result.full, result.drop = result.drop, lod = result.full[1,4]) )

}


