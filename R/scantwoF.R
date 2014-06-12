#' Two-dimensional genome scan with a two-QTL model for function valued trait
#' data.
#'
#' Extension of the R/qtl function \code{\link[qtl]{scantwo}}. Two-dimensional genome
#' scan with a two-QTL model for function valuded trait data.
#' 
#'
#' @param cross An object of class \code{"cross"}. See \code{\link[qtl]{read.cross}} for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param usec Which method to use (\code{"slod"} or \code{"mlod"}) in the two-QTL scan.
#' @param n.perm If specified, a permutation test is performed rather than an
#' analysis of the observed data.  This argument defines the number of
#' permutation replicates.
#' @param \dots More parameters controled in 'scantwo'. See 'scantwo' for
#' details.
#' @export
#' @return If \code{n.perm} is missing, the function returns a list with class
#' \code{"scantwo"} and containing three components.  The first component is a
#' matrix of dimension [tot.pos x tot.pos]; the upper triangle contains the
#' (S/MLOD) scores for the additive model, and the lower triangle contains the
#' LOD scores for the full model.  The diagonal contains the results of
#' a single-QTL scan. The second component of the output is a data.frame indicating the
#' locations at which the two-QTL (S/MLOD) scores were calculated.  The first
#' column is the chromosome identifier, the second column is the position in
#' cM, the third column is a 1/0 indicator for ease in later pulling out only
#' the equally spaced positions, and the fourth column indicates whether the
#' position is on the X chromosome or not.  The final component is a version of
#' the results of 'scanone' including sex and/or cross direction as additive
#' covariates, which is needed for a proper calculation of conditional (S/MLOD)
#' scores.
#'
#' If \code{n.perm} is specified, the function returns a list with six different LOD
#' scores from each of the permutation replicates. ... need more ..
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link[qtl]{scantwo}}, code{\link[funqtl]{scanoneF}}
#' @keywords models
#' @examples
#' data(exd)
#' exd <- calc.genoprob(exd, step = 0)
#' out <-scantwoF(exd, method = "hk", usec="slod")
#' plot(out)
#'
#' # Permutation tests
#' n.perm <- 1000
#' \dontrun{permo <- scantwoF(exd, method="hk", n.perm=n.perm)
#' summary(permo, alpha=0.05)}
scantwoF <- function(cross, pheno.cols, usec=c("slod","mlod"), n.perm, ...) {

    n = nind(cross)
    usec <- match.arg(usec)

    if (missing(pheno.cols))
    pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
    stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    if (missing(n.perm))
    n.perm <- 0

    if (n.perm > 0 ) {
        temp <- cross
        pheno <- as.data.frame(cross$pheno[, pheno.cols])

        Slod <- NULL;
        Mlod <- NULL;
        SlodsH <- NULL;
        SlodsL <- NULL;
        MlodsH <- NULL;
        MlodsL <- NULL;
        Slods <- NULL;
        Mlods <- NULL;

        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]

            out2 <- scantwo(temp, ...)
            out1 <- scanone(temp, ...)

            # out3 for slod
            out3 <- out2
            out3$lod <- apply(out2$lod, 1:2, mean)

            # out4 for mlod
            out4 <- out2
            out4$lod <- apply(out2$lod, 1:2, max)

            SLOD <- rowMeans(out1[,-(1:2)])
            MLOD <- apply(out1[,-(1:2)], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )
            SlodsH <- c(SlodsH, max(summary(out3)$lod.int) )
            SlodsL <- c(SlodsL, max(summary(out3)$lod.fv1) )
            MlodsH <- c(MlodsH, max(summary(out4)$lod.int) )
            MlodsL <- c(MlodsL, max(summary(out4)$lod.fv1) )
        }


#        out <- list(one_slod = as.matrix(Slods), one_mlod = Mlods, fullvadd_slod = SlodsH, fullvadd_mlod = MlodsH, fv1_slod = SlodsL, fv1_mlod = MlodsL)
        out <- list(one = cbind(Slods, Mlods), fullvadd = cbind(SlodsH, MlodsH), fv1 = cbind(SlodsL, MlodsL) )
        class(out) <- c("scantwoperm", "list")
        out

    } else {
        out <- scantwo(cross, pheno.col = pheno.cols, ...)

        if(usec=="slod") {
            out$lod <- apply(out$lod, 1:2, mean)
        }
        if(usec=="mlod") {
            out$lod <- apply(out$lod, 1:2, max)
        }
        out
    }
}
