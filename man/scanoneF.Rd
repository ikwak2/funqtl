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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cross, pheno.cols, n.perm, ...)
{
    n = nind(cross)
    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)
    if (!all(pheno.cols \%in\% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))
    if (missing(n.perm))
        n.perm <- 0
    if (n.perm > 0) {
        cross <- calc.genoprob(cross, step = 0)
        temp <- cross
        pheno <- cross$pheno
        Slods <- NULL
        Mlods <- NULL
        for (rep in 1:n.perm) {
            temp$pheno <- pheno[sample(n), ]
            out <- scanone(temp, pheno.col = pheno.cols, ...)
            SLOD <- rowMeans(out[, -(1:2)])
            MLOD <- apply(out[, -(1:2)], 1, max)
            Slods <- c(Slods, max(SLOD))
            Mlods <- c(Mlods, max(MLOD))
        }
        return(rbind(Slods, Mlods))
    }
    else {
        out <- scanone(cross, pheno.col = pheno.cols, ...)
        SLOD <- rowMeans(out[, -(1:2)])
        MLOD <- apply(out[, -(1:2)], 1, max)
        out[, 3] <- SLOD
        out[, 4] <- MLOD
        names(out)[3:4] <- c("slod", "mlod")
        out[, 1:4]
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
