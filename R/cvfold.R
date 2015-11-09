#' Do cross validation to decide the number of basis function to use in smoothing.
#'
#'
#' @param cross An object of class \code{"cross"}. See the \code{\link[qtl]{read.cross}} for details.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param basisset The set of basis numbers to evaluate.
#' @param fold The number of folder in cross validation.
#' @param random randomly divide folder on times if TRUE and select folders equily spaced time points if FALSE.
#' @return It gives a vector of sum of squared erros for each basis set.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{calcfunpca}}
#' @keywords utilities
#' @export
#' @importFrom fda create.bspline.basis fd pca.fd eval.fd eval.basis
#' @examples
#' data(exd)
#' exd <- calc.genoprob(exd, step=2)
#' cvout <- cvfold(exd, basisset = 4:7, fold = 10)
#' cvout # basis number 5 have the smallest sse. So we take nbasis = 5.
#' Y <- calcfunpca(exd, criteria=0.9, nbasis = 5)$Y
#' out1 <- scanoneM(exd, Y, method = "hk")


cvfold <- function (cross, pheno.cols, basisset, fold = 10, random = TRUE )
{
    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)
    Y = t(cross$pheno[, pheno.cols])
    hasmissing <- apply(Y, 2, function(a) any(is.na(a)))
    if (any(hasmissing))
        Y <- Y[,!hasmissing ]


    m = nrow(Y)

    if(random == TRUE)
    {
        o <- sample(m)
    } else {

        AA <- matrix(1:(floor(m/fold)*fold), nrow = fold)

        o <- NULL
        for (i in 1:fold)
            o <- c(o, AA[i, ])
        o <- c(o, (1:m)[-(1:(floor(m/fold)*fold))] )
    }

    setnum <- floor(m/fold)

    tt <- 0:(m - 1) + 0.5

    tsterr = NULL;
    for(bsis in basisset)
    {
        sumerr <- 0
        for( i in 1:fold)
        {
            foldindex <- o[(1:setnum)+setnum*(i-1)]

            testset <- Y[foldindex,]
            evalset <- Y[-foldindex,]

            testtime <- tt[foldindex]
            evaltime <- tt[-foldindex]

            splinebasis.y <- create.bspline.basis(c(0, m), nbasis = bsis, 4)
            mat <- eval.basis(evaltime, splinebasis.y)
            coef.y <- solve(crossprod(mat), crossprod(mat, evalset))
            yfd = fd(coef.y, splinebasis.y, list("time", "indv", "value"))

            sumerr = sumerr + sum(eval.fd(testtime, yfd) - testset)^2
        }
        tsterr <- c(tsterr, sumerr/fold )
    }
    names(tsterr) <- basisset
    tsterr
}
