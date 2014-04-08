#' QTL scan using dimensional reduction.
#'
#' Genome scan with a single QTL model using dimensional reduction
#'
#' @param cross An object of class \code{"cross"}. See \code{\link[qtl]{read.cross}} for details.
#' @param Y Dimension-reduced data set.
#' @param tol Tolerance; passed to \code{\link[stats]{lm.fit}}
#' @param n.perm If specified, a permutation test is performed rather than an
#' analysis of the observed data.  This argument defines the number of
#' permutation replicates.
#' @param method The \code{"hk"} option use multi-trait QTL mapping proposed by Haley
#' and Knott. The \code{"f"} option use an FLOD score.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @return If \code{n.perm} is missing, the function returns a data.frame whose
#' first two columns contain the chromosome IDs and cM positions.  Subsequent
#' third and fourth columns contain the SLOD and MLOD scores.
#'
#' If \code{n.perm} is specified, the function returns the results of a permutation
#' test and the output returns the matrix of two columns. The first column for
#' SLOD and the second column for MLOD score.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link[qtl]{scanone}}, \code{\link{scanoneF}}
#' @keywords model
#' @export
#' @examples
#' data(simspal)
#'
#' \dontshow{simspal <- subset(simspal,chr=1:4,ind=1:100)}
#' # Genotype probabilities for H-K
#' simspal <- calc.genoprob(simspal)
#'
#' # dimensional reduction of Y
#' Y <- calcpca(simspal, criteria=.999)
#' 
#' # do multitrait mapping
#' out.hk <- scanoneM(simspal, Y=Y, method="hk")
#' out.f  <- scanoneM(simspal, Y=Y, method="f")
#' out.sl <- scanoneM(simspal, Y=Y, method="sl")
#' out.ml <- scanoneM(simspal, Y=Y, method="ml")
#'
#' # Summarize results
#' summary(out.hk)
#' summary(out.f)
#' summary(out.sl)
#' summary(out.ml)
#' 
#' # Plot the results
#' par(mfrow=c(3,1))
#' plot(out.hk)
#' plot(out.f)
#' plot(out.sl, out.ml)

scanoneM <- function(cross, Y, tol=1e-7, n.perm=0, method=c("hk","f", "sl", "ml"), pheno.cols ) {

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)
    n.ind <- nind(cross)
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

            if(method == "sl") {
                return(outs[,1:3])
            } else {
                return(outs[,c(1,2,4)])
            }

        } else {

            E <- matrix(NA, n.ind, p)
            X <- cbind(rep(1,n.ind))
            E <- lm.fit(X, Y, tol=tol)$residuals
            Sigma <- crossprod(E)

            if( method == "hk") {
                L0 <- determinant(Sigma)$modulus
            } else {
                L0 <- sum(diag(Sigma))
            }

            out <- NULL;
            for(i in 1:n.chr) {
                LOD = NULL;
                map <- attr(cross$geno[[i]]$prob, "map")
                for(j in 1:length(map)) {
                    X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                    E <- lm.fit(X, Y, tol=tol)$residuals
                    Sigma <- crossprod(E)

                    if( method == "hk") {
                        L1 <- determinant(Sigma)$modulus
                    } else {
                        L1 <- sum(diag(Sigma))
                    }

                    LOD <- c(LOD, n.ind/2*log(L0/L1,10 ) )
                }
                out <- rbind(out, cbind(rep(as.numeric(chrnames(cross)[i]),length(map)), map, LOD) )
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
                cat("Permutation", rep,"\n")
                o <- sample(n.ind)
                if(is.vector(Y)) { nY <- Y[o] } else { nY <- Y[o,] }


                E <- matrix(NA, n.ind, p)
                X <- cbind(rep(1,n.ind))
                E <- lm.fit(X, nY, tol=tol)$residuals
                Sigma <- crossprod(E)

                if( method == "hk") {
                    L0 <- determinant(Sigma)$modulus
                } else {
                    L0 <- sum(diag(Sigma))
                }

                out <- NULL;
                for(i in 1:n.chr) {
                    LOD = NULL;
                    map <- attr(cross$geno[[i]]$prob, "map")
                    for(j in 1:length(map)) {
                        X <- cbind(rep(1,n.ind), cross$geno[[i]]$prob[,j,1])
                        E <- lm.fit(X, nY, tol=tol)$residuals
                        Sigma <- crossprod(E)

                        if( method == "hk") {
                            L1 <- determinant(Sigma)$modulus
                        } else {
                            L1 <- sum(diag(Sigma))
                        }
                        LOD <- c(LOD, n.ind/2*log(L0/L1,10 ) )

                    }
                    out <- rbind(out, cbind(rep(as.numeric(chrnames(cross)[i]),length(map)), map, LOD) )
                }

                lods <- c(lods, max(out[,3]) )
            }
            class(lods) <- c("scanoneperm","vector")
            return(lods)
        }
    }
}



