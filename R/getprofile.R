#' Make a profile LOD score matrix for a given multiple QTL model.
#'
#' Calculate profile LOD scores for each chromosome.
#'
#'
#' @param cross
#'
#' An object of class 'cross'. See 'read.cross' for detail.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param qtl A QTL object, as produced by 'makeqtl', containing the positions
#' of the QTL.  Provide either 'qtl' or the pair 'chr' and 'pos'.
#' @param chr Vector indicating the chromosome for each QTL; if 'qtl' is
#' provided, this should not be.
#' @param pos Vector indicating the positions for each QTL; if 'qtl' is
#' provided, this should not be.
#' @param qtl.name Optional user-specified name for each QTL.  If 'qtl' is
#' provided, this should not be.
#' @param covar A matrix or data.frame of covariates.  These must be strictly
#' numeric.
#' @param formula An object of class 'formula' indicating the model to be
#' fitted.  (It can also be the character string representation of a formula.)
#' QTLs are indicated as 'Q1', 'Q2', etc.  Covariates are indicated by their
#' names in 'covar'.
#' @param method Indicates whether to use multiple imputation or Haley-Knott
#' regression.
#' @param model The phenotype model: the usual model or a model for binary
#' traits
#' @param verbose If TRUE, give feedback about progress.  If 'verbose' is an
#' integer > 1, further messages from 'scanqtl' are also displayed.
#' @param tol Tolerance for convergence for the binary trait model.
#' @param maxit.fitqtl Maximum number of iterations for fitting the binary
#' trait model.
#' @param tpy type of output. If there are more QTL's in one chromosome. It
#' plot them separately if tpy = "sep", On the other hand, it combine then in
#' one chromosome taking maximum values of them if tpy = "comb".
#' @return A 'lodprofileM' or 'lodprofileM2' object.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{plotprofile}}
#' @keywords models
#' @examples
#'
#' data(simspal)
#' simspal <- calc.genoprob(simspal, step=0)
#'
#' # difference between tpy = "sep" and "comb"
#'
#' par(mfrow=c(1,2))
#' qtlslod <- makeqtl(simspal, chr = c(1, 1, 4),
#'                pos = c(36.6, 61, 27.8), what = "prob")
#' lodmat1 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =1:241,
#'                     formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "sep")
#' lodmat2 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =1:241,
#'                     formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "comb")
#'
#' plotprofile(lodmat1, main="The Profile LOD image of data")
#' plotprofile(lodmat2, main="The Profile LOD image of data")
#'

getprofile <-
function (tpy=c("comb","sep"), ...) {

    tpy <- match.arg(tpy)

    if( tpy == "comb" ) {
        out <- profileLodMatfn2(...)
    } else if (tpy == "sep") {
        out <- profileLodMatfn(...)
    }
    out
}
