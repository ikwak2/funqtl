#' Do dimensional reduction using pca.
#'
#' Do dimensional reduction using pca.
#'
#'
#' @param cross An object of class 'cross'. See 'read.cross' for details.
#' @param n.max The number of maximum reduced dimension.
#' @param criteria how much of variance explained.
#' @param nn The number of exact reduced dimension
#' @return It gives a matrix that each column have principal components.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{scanoneM}}
#' @keywords utilities
#' @examples
#'
#'      data(exd)
#'
#'      exd <- calc.genoprob(exd, step=2)
#'      Y <- calcpca(exd, criteria=.9)
#'      out1 <- scanoneM(exd, Y, method = "hk")
#'      out2 <- scanoneM(exd, Y, method = "f")
#'      out3 <- scanoneM(exd, Y, method = "sl")
#'      out4 <- scanoneM(exd, Y, method = "ml")
#'

calcpca <- function(cross, n.max=5, criteria=.9, nn = 0) {

    Y = cross$pheno
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
    pc <- udv$u %*% diag(udv$d)
    pc[,1:nn]
}