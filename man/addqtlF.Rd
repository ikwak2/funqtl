\name{addqtlF}
\alias{addqtlF}

\title{
Search for Multiple QTL model
}
\description{
Scan for an additional QTL for function valued trait data set. Modified
version of addqtl function in qtl package.
}
\usage{
addqtlF(cross, pheno.cols, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
}
  \item{\dots}{
    See 'addqtl' for details.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

An object of class 'scanone', as produced by the 'scanone'
function.  LOD scores are relative to the base model (with any
terms that include the new QTL omitted).
%%                        ##    object have slod and mlod scores in 3rd
%%                        and 4th column is this ok?? ###

}
\references{
     Haley, C. S. and Knott, S. A. (1992) A simple regression method
     for mapping quantitative trait loci in line crosses using flanking
     markers.  _Heredity_ *69*, 315-324.

     Sen, S. and Churchill, G. A. (2001) A statistical framework for
     quantitative trait mapping.  _Genetics_ *159*, 371-387.

}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
     'addqtl', 'scanone', 'scanoneF', 'fitqtl', 'scanqtl', 'refineqtl', 'refineqtlF',
     'makeqtl', 'addtoqtl', 'addpair', 'addint'
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cross, pheno.cols, ...)
{
    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)
    if (!all(pheno.cols \%in\% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))
    LODS <- NULL
    for (i in pheno.cols) {
        out <- addqtl(cross, pheno.col = i, ...)
        LODS <- cbind(LODS, out$lod)
    }
    MXy <- max(LODS)
    Slods <- apply(LODS, 1, mean)
    Mlods <- apply(LODS, 1, max)
    out[, 3] <- Slods
    out[, 4] <- Mlods
    names(out)[3:4] <- c("slod", "mlod")
    out[, 1:4]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
