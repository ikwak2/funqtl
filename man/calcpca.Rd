\name{calcpca}
\alias{calcpca}
\title{Do dimensional reduction using pca.}
\usage{
calcpca(cross, pheno.cols, n.max = 5, criteria = 0.9, nn = 0)
}
\arguments{
  \item{cross}{An object of class \code{"cross"}. See the
  \code{\link[qtl]{read.cross}} for details.}

  \item{pheno.cols}{Columns in the phenotype matrix to be
  used as the phenotype.}

  \item{n.max}{The number of maximum reduced dimension.}

  \item{criteria}{how much of variance explained.}

  \item{nn}{The number of exact reduced dimension}
}
\value{
It gives a matrix that each column have principal
components.
}
\description{
Do dimensional reduction using pca.
}
\examples{
data(exd)
exd <- calc.genoprob(exd, step=2)
Y <- calcpca(exd, criteria=0.9)
out1 <- scanoneM(exd, Y, method = "hk")
}
\author{
Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}
\seealso{
\code{\link{scanoneM}}
}
\keyword{utilities}

