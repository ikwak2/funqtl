#' Estimate QTL effects
#'
#' Estimate QTL effect at each genomic position for each of multiple phenotypes.
#' However, this function uses a very simple method to predict the effect, and it is recommended to use the fitqtl function referring to a [vignette](https://ikwak2.github.io/tutorials/funqtl.html) to actually predict the effect.
#'
#' @param cross An object of class 'cross'. See 'read.cross' for details.
#' @param pheno.cols phenotype columns to be used
#' @return A matrix of coefficients. (i,j)th item is a coefficient of jth
#' position as a qtl of ith observation.
#' @author Karl W Broman, Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{plotlod}}
#' @keywords models
#' @export
#' @importFrom stats lm
#' @examples
#' data(simspal)
#' simspal <- calc.genoprob(simspal)
#' phe <- 1:nphe(simspal)
#' \dontshow{phe <- seq(1, nphe(simspal), by=60)}
#' eff <- geteffects(simspal, pheno.cols=phe)
#'
geteffects <- function(cross,pheno.cols) {
    if(missing(pheno.cols))
        pheno.cols=1:nphe(cross)
    if( len(attr(cross$geno[[1]]$prob, "dimnames")[[3]]) != 2 )
        stop("This function works for backcross or recombinant inbred lines.")

    
    phe <- as.matrix(cross$pheno[,pheno.cols])
    eff <- NULL
    for(i in 1:nchr(cross)) {
        pr <- cross$geno[[i]]$prob[,,2]
        eff <- rbind(eff, t(apply(pr, 2, function(a,b) lm(b~a)$coef[2,], phe)))
    }
    eff
}
