#' Genome scan with a single QTL model
#'
#' Extension of 'scanone' function of 'qtl' package. Genome scan with a single
#' QTL model, with possible allowance for covariates, using several possible
#' models for the function valued phenotype.
#'
#'
#' @param cross An object of class 'cross'. See 'read.cross' for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param n.perm If specified, a permutation test is performed rather than an
#' analysis of the observed data.  This argument defines the number of
#' permutation replicates.
#' @param \dots More parameters controled in 'scanone'. See 'scanone' for
#' details.
#' @return If 'n.perm' is missing, the function returns a data.frame whose
#' first two columns contain the chromosome IDs and cM positions.  Subsequent
#' third and fourth columns contain the SLOD and MLOD scores.
#'
#' If 'n.perm' is specified, the function returns the results of a permutation
#' test and the output returns the matrix of two columns. The first column for
#' SLOD and the second column for MLOD score.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
#' \code{\link[qtl]{summary.scanone}}, \code{\link[qtl]{calc.genoprob}}
#' @keywords models
#' @examples
#'
#'
#'      data(exd)
#'
#'      ## genotype probabilities
#'      exd <- calc.genoprob(exd, step=2)
#'      out <- scanoneF(exd)
#'
#'      # Summarize results
#'      summary(out)
#'
#'      # Plot the results : red for slod, blue for mlod
#'      plot(out, out[,c(1,2,4)], col = c("red","blue"), ylab = "lod")
#'
#'
#'      # Permutation tests
#' \dontrun{perm1 <- scanoneF(exd, method="hk", n.perm=1000)}
#' \dontshow{perm1 <- scanoneF(exd, method="hk", n.perm=5)}
#'      summary(perm1, alpha=c(0.05, 0.10))
#'
#'      # Results above the 0.05 threshold
#'      summary(out, perms=perm1, alpha=0.05)
#'
#'

scanoneF <-
function(cross, pheno.cols, n.perm, ...) {

    n = nind(cross)

    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    if (missing(n.perm))
        n.perm <- 0

    if (n.perm > 0 ) {

        temp <- cross
        pheno <- cross$pheno[,pheno.cols]

        Slods <- NULL;
        Mlods <- NULL;
        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]

            out <- scanone(temp, pheno.col = pheno.cols, ...)
            SLOD <- rowMeans(out[,-(1:2),drop=FALSE])
            MLOD <- apply(out[,-(1:2),drop=FALSE], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )

        }

        permout <- cbind(Slods,Mlods)
        colnames(permout) <- c("slod","mlod")
        class(permout) <- c("scanoneperm","matrkx")


        return( permout )

    } else {

        out <- scanone(cross, pheno.col = pheno.cols, ...)
        SLOD <- rowMeans(out[,-(1:2),drop=FALSE])
        MLOD <- apply(out[,-(1:2),drop=FALSE], 1, max)

        out[,3] <- SLOD
        out[,4] <- MLOD
        names(out)[3:4] <- c("slod","mlod")
        out[,1:4]
    }
}
