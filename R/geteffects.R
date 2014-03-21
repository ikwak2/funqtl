#' get coefficients of qtl. (get sign information)
#'
#' required function to run plotlod which plot signed lod image. The plotlod
#' function need a sign information for each LOD score. This geteffects
#' function get coefficients at each gene position. So we can get sign
#' information from this function.
#'
#'
#' @param cross An object of class 'cross'. See 'read.cross' for details.
#' @param pheno.cols phenotype columns to be used
#' @return A matrix of coefficients. (i,j)th item is a coefficient of jth
#' position as a qtl of ith observation.
#' @author Karl W Broman, Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso 'plotlod'
#' @keywords models
#' @examples
#'
#'
#' data(simspal)
#' simspal <- calc.genoprob(simspal, step=1)
#'
#' out <- scanone(simspal, pheno.col=1:241 , method="hk")
#' eff <- geteffects(simspal, pheno.cols=1:241)
#'
#' plotlod(out, eff, gap=15)
#'

geteffects <- function(cross,pheno.cols) {
    if(missing(pheno.cols))
        pheno.cols=1:nphe(cross)

#    out <- scanone(cross, pheno.col = pheno.cols, method="hk")
    phe <- as.matrix(cross$pheno)
    eff <- NULL
    for(i in 1:nchr(cross)) {
        pr <- cross$geno[[i]]$prob[,,2]
        eff <- rbind(eff, t(apply(pr, 2, function(a,b) lm(b~a)$coef[2,], phe)))
    }
    eff
}
