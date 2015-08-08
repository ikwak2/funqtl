\name{cvfold}
\alias{cvfold}
\title{Do cross validation to decide the number of basis function to use in smoothing.}
\usage{
cvfold(cross, pheno.cols, basisset, fold = 10, random = TRUE)
}
\arguments{
  \item{cross}{An object of class \code{"cross"}. See the
  \code{\link[qtl]{read.cross}} for details.}

  \item{pheno.cols}{Columns in the phenotype matrix to be
  used as the phenotype.}

  \item{basisset}{The set of basis numbers to evaluate.}

  \item{fold}{The number of folder in cross validation.}

  \item{random}{randomly divide folder on times if TRUE and select folders equily spaced time points if FALSE.}
}
\value{
It gives a vector of sum of squared erros for each basis
set.
}
\description{
Do cross validation to decide the number of basis function
to use in smoothing.
}
\examples{
data(exd)
exd <- calc.genoprob(exd, step=2)
cvout <- cvfold(exd, basisset = 4:7, fold = 10)
cvout # basis number 5 have the smallest sse. So we take nbasis = 5.
Y <- calcfunpca(exd, criteria=0.9, nbasis = 5)
out1 <- scanoneM(exd, Y, method = "hk")
}
\author{
Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}
\seealso{
\code{\link{calcfunpca}}
}
\keyword{utilities}
