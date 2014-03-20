\name{scantwoF}
\alias{scantwoF}
\title{
  Two-dimensional genome scan with a two-QTL model for function valued
  trait data.
}
\description{
Extension of 'scantwo' function of 'qtl' package.
Two-dimensional genome scan with a two-QTL model for function valuded
trait data.
}
\usage{
scantwoF(cross, pheno.cols, usec=c("slod","mlod"), n.perm, ...)
}
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
}
  \item{usec}{
    Which method to use between 'slod' and 'mlod' in two-QTL scan.
}
  \item{n.perm}{
    If specified, a permutation test is performed rather than an
    analysis of the observed data.  This argument defines the
    number of permutation replicates.

}
  \item{\dots}{
    More parameters controled in 'scantwo'. See 'scantwo' for details.
}
}
\value{
  If 'n.perm' is missing, the function returns a list with class
     '"scantwo"' and containing three components.  The first component
     is a matrix of dimension [tot.pos x tot.pos]; the upper triangle
     contains the (S/MLOD) scores for the additive model, and the lower
     triangle contains the LOD scores for the full model.  The diagonal
     contains the results of 'scanone'. The second component of the
     output is a data.frame indicating the locations at which the
     two-QTL (S/MLOD) scores were calculated.  The first column is the
     chromosome identifier, the second column is the position in cM,
     the third column is a 1/0 indicator for ease in later pulling out
     only the equally spaced positions, and the fourth column indicates
     whether the position is on the X chromosome or not.  The final
     component is a version of the results of 'scanone' including sex
     and/or cross direction as additive covariates, which is needed for
     a proper calculation of conditional (S/MLOD) scores.

     If 'n.perm' is specified, the function returns a list with six
     different LOD scores from each of the permutation replicates.
... need more ..
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

\seealso{
  \code{\link[qtl]{scantwo}}, \code{\link[qtl]{summary.scantwo}}, \code{\link[qtl]{calc.genoprob}}, \code{\link[funqtl]{scanoneF}}
}
\examples{

data(exd)
exd <- calc.genoprob(exd, step = 0)

out <-scantwoF(exd, method = "hk", usec="slod")
plot(out)

# permutations

# Permutation tests
\dontrun{permo <- scantowF(exd, method="hk", n.perm=1000)}
\dontshow{permo <- scanone(exd, method="hk", n.perm=5)}

summary(permo, alpha=0.05)




}
\keyword{ models }
