#' Do dimensional reduction using pca.
#'
#' Do dimensional reduction using pca.
#'
#'
#' @param cross An object of class \code{"cross"}. See the \code{\link[qtl]{read.cross}} for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param n.max The number of maximum reduced dimension.
#' @param criteria how much of variance explained.
#' @param nn The number of exact reduced dimension
#' @return It gives a matrix that each column have principal components.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{scanoneM}}
#' @keywords utilities
#' @export
#' @examples
#' data(exd)
#' exd <- calc.genoprob(exd, step=2)
#' Y <- calcpca(exd, criteria=0.9)
#' out1 <- scanoneM(exd, Y, method = "hk")
calcpca <- function(cross, pheno.cols, n.max=5, criteria=.9, nn = 0) {

    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    Y = cross$pheno[, pheno.cols]
    udv <- svd(Y)
    vec <- udv$d^2
    vec <- vec / sum(vec)

    ss = 0

    if( nn == 0) {
        for(j in 1:n.max) {
            ss = ss + vec[j]
            if ( ss > criteria ) {
                nn = j
                break
            }
        }
    }

    if (j == n.max & ss < criteria) {
        stop("You should incarese n.max to meet your criteria.")
    }

    pc <- udv$u %*% diag(udv$d)
    pc[,1:nn]
}
