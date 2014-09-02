#' Do dimensional reduction using functional pca.
#'
#'
#' @param cross An object of class \code{"cross"}. See the \code{\link[qtl]{read.cross}} for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param n.max The number of maximum reduced dimension.
#' @param criteria how much of variance explained.
#' @param nbasis The number of basis to use.
#' @param nn The number of exact reduced dimension
#' @return It gives a matrix that each column have principal components.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{scanoneM}}
#' @keywords utilities
#' @export
#' @importFrom fda create.bspline.basis fd pca.fd eval.fd
#' @examples
#' data(exd)
#' exd <- calc.genoprob(exd, step=2)
#' cvout <- cvfold(exd, basisset = 4:7, fold = 10)
#' cvout # basis number 5 have the smallest sse. So we take nbasis = 5.
#' Y <- calcfunpca(exd, criteria=0.9, nbasis = 5)
#' out1 <- scanoneM(exd, Y, method = "hk")

calcfunpca <- function(cross, pheno.cols, n.max=4, criteria=.9, nbasis, nn = 0) {

    if(missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if(n.max > nbasis) 
        n.max = nbasis

    Y = t( cross$pheno[,pheno.cols] )

    hasmissing <- apply(Y, 1, function(a) any(is.na(a)) )
    if(any(hasmissing) )
        Y <- Y[!hasmissing, ]

    m = nrow(Y)

    if(missing(nbasis) )
       nbasis = m
    splinebasis.y <- create.bspline.basis(c(0,m), nbasis, 4)

    time <- 0:(m-1) + 0.5
    mat <- eval.basis(time, splinebasis.y)

    coef.y <- solve(crossprod(mat), crossprod(mat, Y))

    yfd = fd(coef.y, splinebasis.y, list("time","indv","value") )

    if( nn == 0) {
        y.pcalist3 = pca.fd(yfd, n.max)

        ss = 0
        for(j in 1:n.max) {
            ss = ss + y.pcalist3$varprop[j]
            if ( ss > criteria ) {
                nn = j
                break
            }
        }
        if( j == n.max & ss < criteria ) {
            stop("You should increase n.max to meet your criteria.")
        }
	if (j == n.max)
	   nn = n.max
    }

    y.pcalist3 = pca.fd(yfd, nn)
    eigfc3 <- y.pcalist3$harmonics
    mat3 <- eval.fd(time, eigfc3)
    nY3 <- t(solve(crossprod(mat3), crossprod(mat3, Y) ) )
}
