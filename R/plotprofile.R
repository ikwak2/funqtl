#' Plot profile LOD maxrix.
#'
#' Plot profile LOD maxrix of given multiple QTL model.
#'
#'
#' @param lodmatlist An 'lodprofileM' object that produced by 'profileLodMatfn'
#' function.
#' @param ylab A lable of y. Default is "QTL position".
#' @param xlab A lable of x. Default is "Time".
#' @param mval The maximum LOD value of legend. The color of legend goes 0 to
#' 'mval'. If this value is less than the max lod score, it automatically
#' changed to max value.
#' @param \dots More graphical components of 'image.plot'.
#' @return A graph of profile LOD matrix.
#' @author Il-Youp Kwak, <ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{getprofile}}
#' @keywords hplot
#' @examples
#'
#' data(simspal)
#' simspal <- calc.genoprob(simspal, step=0)
#'
#' # difference between tpy = "sep" and "comb"
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

plotprofile <- function(lodmatlist, ylab="QTL position", xlab="Time", mval=0, ...) {

    if(class(lodmatlist) == "lodprofileM") {
        plotlodmatlist(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    } else if (class(lodmatlist) == "lodprofileM2") {
        plotlodmatlist2(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    }
}
