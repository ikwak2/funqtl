\name{scanoneM}
\alias{scanoneM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Genome scan with a single QTL model using dimensional reduction
}
\description{
  QTL scan using dimensional reduction. 
}
\usage{
scanoneM(cross, Y, pheno.cols, n.perm, method, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
}
  \item{Y}{
    Demension reduced data set. getY(cross) get reduced data set using PCA.
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
     Y <- getY(exd, criteria=.9)
     out1 <- scanoneM(exd, Y, method = "hk")
     out2 <- scanoneM(exd, Y, method = "f")
     out3 <- scanoneM(exd, Y, method = "sl")
     out4 <- scanoneM(exd, Y, method = "ml")
 

     # Summarize results
     summary(out)
     
     # Plot the results : red for slod, blue for mlod
     plot(out, out[,c(1,2,4)], col = c("red","blue"), ylab = "lod")
     

     # Permutation tests
     perm1 <- scanone(exd, method="hk", n.perm=100)
     # Threshold from the permutation test
     summary(perm1, alpha=c(0.05, 0.10))
     
     # Results above the 0.05 threshold
     summary(out, perms=perm1, alpha=0.05)
     
     
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
