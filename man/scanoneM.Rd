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
  \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
  \code{\link[qtl]{summary.scanone}},   \code{\link[qtl]{calc.genoprob}}
}
\examples{

data(simspal)

\dontshow{simspal <- subset(simspal,chr=1:4,ind=1:100)}
# Genotype probabilities for H-K
\dontrun{simspal <- calc.genoprob(simspal, step=2.5)
}\dontshow{simspal <- calc.genoprob(simspal)}

# dimensional reduction of Y
Y <- calcpca(simspal, criteria=.999)

# do multitrait mapping
out.hk <- scanoneM(simspal, Y=Y, method="hk")
out.f  <- scanoneM(simspal, Y=Y, method="f")
out.sl <- scanoneM(simspal, Y=Y, method="sl")
out.ml <- scanoneM(simspal, Y=Y, method="ml")

# Summarize results
summary(out.hk)
summary(out.f)
summary(out.sl)
summary(out.ml)


# Plot the results
par(mfrow=c(3,1))
plot(out.hk)
plot(out.f)
plot(out.sl, out.ml)


}

\keyword{ model }
