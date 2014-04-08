#' Search for Multiple QTL model
#'
#' Scan for an additional QTL for function valued trait data set. Modified
#' version of addqtl function in qtl package.
#'
#'
#' @param cross An object of class \code{"cross"}. See \code{\link[qtl]{read.cross}} for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param \dots Passed to the R/qtl function \code{\link[qtl]{addqtl}}.
#' @return An object of class \code{"scanone"}, as produced by the R/qtl function \code{\link[qtl]{scanone}}.
#' LOD scores are relative to the base model (with any terms that include the
#' new QTL omitted).
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link[qtl]{addqtl}}, \code{\link[qtl]{scanone}},
#' \code{\link{scanoneF}}, \code{\link[qtl]{makeqtl}}
#' @references Haley, C. S. and Knott, S. A. (1992) A simple regression method
#' for mapping quantitative trait loci in line crosses using flanking markers.
#' _Heredity_ *69*, 315-324.
#' @keywords models
#' @export
#' @examples
#'
#' data(simspal)
#' \dontshow{simspal <- subset(simspal, chr=c(1,4), ind=1:50)}
#'
#' qtl1 <- makeqtl(simspal, chr = 4, 27.8, what="prob")
#'
#' phe <- 1:nphe(simspal)
#' \dontshow{phe <- seq(1, nphe(simspal), by=60)}
#' added <- addqtlF(simspal, qtl =  qtl1, pheno.cols =phe, method="hk")
#' max(added)
#' plot(added)
addqtlF <-
function(cross, pheno.cols, ...) {

    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    LODS <- NULL;
    for(i in pheno.cols ) {
        out <- addqtl(cross, pheno.col = i, ...)
        LODS <- cbind(LODS, out$lod)
    }

    MXy <- max(LODS)
    Slods <- apply(LODS, 1, mean)
    Mlods <- apply(LODS, 1, max)
    out[,3] <- Slods
    out[,4] <- Mlods
    names(out)[3:4] <- c("slod","mlod")

    out[,1:4]
}
