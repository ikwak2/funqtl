\name{scanoneM}
\alias{scanoneM}

\title{
  Genome scan with a single QTL model using dimensional reduction
}
\description{
  QTL scan using dimensional reduction.
}
\usage{
scanoneM(cross, Y, tol=1e-7, n.perm=0, method=c("hk","f", "sl", "ml"), pheno.cols )

}

\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{Y}{
    Demension reduced data set. getY(cross) get reduced data set using PCA.
}
\item{tol}{
  Tolerance; passed to lm.fit
}
  \item{n.perm}{
    If specified, a permutation test is performed rather than an
    analysis of the observed data.  This argument defines the
    number of permutation replicates.
}

\item{method}{
  The "hk" option use multi-trait QTL mapping propsoed by Haley and
  Knott. The "f" option use FLOD score.
}
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
}

}
\value{
     If 'n.perm' is missing, the function returns a data.frame whose
     first two columns contain the chromosome IDs and cM positions.
     Subsequent third and fourth columns contain the SLOD and MLOD scores.

     If 'n.perm' is specified, the function returns the results of a
     permutation test and the output returns the matrix of two
     columns. The first column for SLOD and the second column for MLOD score.

}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}


\seealso{
  'scanone', 'plot.scanone', 'summary.scanone', 'calc.genoprob',
}
\examples{
     data(exd)

     exd <- calc.genoprob(exd, step=2)
     Y <- calcpca(exd, criteria=.9)
     out1 <- scanoneM(exd, Y, method = "hk")
     out2 <- scanoneM(exd, Y, method = "f")
     out3 <- scanoneM(exd, Y, method = "sl")
     out4 <- scanoneM(exd, Y, method = "ml")
}

\keyword{ model }
