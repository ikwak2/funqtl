\name{scanoneF}
\alias{scanoneF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Genome scan with a single QTL model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Extension of 'scanone' function of 'qtl' package.
Genome scan with a single QTL model, with possible allowance for
covariates, using several possible models for the function valued phenotype.
}
\usage{
scanoneF(cross, pheno.cols, n.perm, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
}
  \item{n.perm}{
    If specified, a permutation test is performed rather than an
    analysis of the observed data.  This argument defines the
    number of permutation replicates.

}
  \item{\dots}{
    More parameters controled in 'scanone'. See 'scanone' for details.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  'scanone', 'plot.scanone', 'summary.scanone', 'calc.genoprob',
}
\examples{

     data(exd)

     exd <- calc.genoprob(exd, step=2)
     out <- scanoneF(exd)

     # Summarize results
     summary(out)

     # Plot the results : red for slod, blue for mlod
     plot(out, out[,c(1,2,4)], col = c("red","blue"), ylab = "lod")


     # Permutation tests
\dontrun{perm1 <- scanoneF(exd, method="hk", n.perm=1000)}
\dontshow{perm1 <- scanoneF(exd, method="hk", n.perm=5)}
     summary(perm1, alpha=c(0.05, 0.10))

     # Results above the 0.05 threshold
     summary(out, perms=perm1, alpha=0.05)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
