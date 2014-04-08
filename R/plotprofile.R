#' Plot profile LOD matrix.
#'
#' Plot profile LOD matrix of given multiple QTL model.
#'
#'
#' @param lodmatlist An \code{"lodprofileM"} object as produced by \code{\link{getprofile}}
#' @param ylab A label of y. Default is "QTL position".
#' @param xlab A label of x. Default is "Time".
#' @param mval The maximum LOD value of legend. The color of legend goes 0 to
#' \code{mval}. If this value is less than the maximum LOD score, it is automatically
#' changed to that maximum.
#' @param \dots Additional graphical components, passed to \code{\link[fields]{image.plot}}.
#' @return None.
#' @author Il-Youp Kwak, <ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{getprofile}}
#' @keywords hplot
#' @export
#' @examples
#' data(simspal)
#' simspal <- calc.genoprob(simspal, step=0)
#' \dontshow{simspal <- subset(simspal,chr=c(1,3,4), ind=1:50)}
#'
#' phe <- 1:nphe(simspal)
#' \dontshow{phe <- seq(1, nphe(simspal), by=60)}
#'
#' # difference between tpy = "sep" and "comb"
#'
#' par(mfrow=c(1,2))
#' qtlslod <- makeqtl(simspal, chr = c(1, 1, 4),
#'                pos = c(36.6, 61, 27.8), what = "prob")
#' lodmat1 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =phe,
#'                     formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "sep")
#' plotprofile(lodmat1, main="tpy=\"sep\"")
#'
#' lodmat2 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =phe,
#'                     formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "comb")
#' plotprofile(lodmat2, main="tpy=\"comb\"")
#'
plotprofile <- function(lodmatlist, ylab="QTL position", xlab="Time", mval=0, ...) {

    if(class(lodmatlist)[1] == "lodprofileM") {
        plotlodmatlist(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    } else if (class(lodmatlist)[1] == "lodprofileM2") {
        plotlodmatlist2(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    }
}
