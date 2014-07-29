\name{scanoneM}
\alias{scanoneM}
\title{QTL scan using dimensional reduction.}
\usage{
scanoneM(cross, Y, tol = 1e-07, n.perm = 0, method = c("hk", "f", "sl",
  "ml"), pheno.cols)
}
\arguments{
  \item{cross}{An object of class \code{"cross"}. See
  \code{\link[qtl]{read.cross}} for details.}

  \item{Y}{Dimension-reduced data set.}

  \item{tol}{Tolerance; passed to
  \code{\link[stats]{lm.fit}}}

  \item{n.perm}{If specified, a permutation test is
  performed rather than an analysis of the observed data.
  This argument defines the number of permutation
  replicates.}

  \item{method}{The \code{"hk"} option use multi-trait QTL
  mapping proposed by Haley and Knott. The \code{"f"}
  option use an FLOD score.}

  \item{pheno.cols}{Columns in the phenotype matrix to be
  used as the phenotype.}
}
\value{
If \code{n.perm} is missing, the function returns a
data.frame whose first two columns contain the chromosome
IDs and cM positions.  Subsequent third and fourth columns
contain the SLOD and MLOD scores.

If \code{n.perm} is specified, the function returns the
results of a permutation test and the output returns the
matrix of two columns. The first column for SLOD and the
second column for MLOD score.
}
\description{
Genome scan with a single QTL model using dimensional
reduction
}
\examples{
data(simspal)

\dontshow{simspal <- subset(simspal,chr=1:4,ind=1:100)}
# Genotype probabilities for H-K
simspal <- calc.genoprob(simspal)

# dimensional reduction of Y
cvout <- cvfold(simspal, basisset = 5:50, fold = 20)
plot(cvout) ## take nbasis = 15
Y <- calcfunpca(simspal, criteria=.999, nbasis = 15)

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
\author{
Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}
\seealso{
\code{\link[qtl]{scanone}}, \code{\link{scanoneF}}
}
\keyword{model}

