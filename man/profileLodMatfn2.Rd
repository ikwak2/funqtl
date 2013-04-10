\name{profileLodMatfn2}
\alias{profileLodMatfn2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
   Make a profile LOD score matrix for a given multiple QTL model.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
   Calculate profile LOD scores for each chromosome. If there are more
   QTL's in one chromosome. It plot them separately, On the other hand, the
   'profileLodMatfn2' combine then in one chromosome taking maximum
   values of them.
}
\usage{
profileLodMatfn2(cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL, formula, method = c("imp", "hk"), model = c("normal", "binary"), verbose = TRUE, tol = 1e-04, maxit.fitqtl = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
An object of class 'cross'. See 'read.cross' for detail.
}
  \item{pheno.cols}{
Columns in the phenotype matrix to be used as the phenotype.
}
  \item{qtl}{
    A QTL object, as produced by 'makeqtl', containing the
    positions of the QTL.  Provide either 'qtl' or the pair 'chr'
    and 'pos'.
}
  \item{chr}{
    Vector indicating the chromosome for each QTL; if 'qtl' is
    provided, this should not be.
}
  \item{pos}{
    Vector indicating the positions for each QTL; if 'qtl' is
    provided, this should not be.
}
  \item{qtl.name}{
    Optional user-specified name for each QTL.  If 'qtl' is
    provided, this should not be.
}
  \item{covar}{
    A matrix or data.frame of covariates.  These must be strictly
    numeric.
}
  \item{formula}{
    An object of class 'formula' indicating the model to be
    fitted.  (It can also be the character string representation
    of a formula.)  QTLs are indicated as 'Q1', 'Q2', etc.
    Covariates are indicated by their names in 'covar'.
}
  \item{method}{
    Indicates whether to use multiple imputation or Haley-Knott
    regression.
}
  \item{model}{
    The phenotype model: the usual model or a model for binary
    traits
}
  \item{verbose}{
    If TRUE, give feedback about progress.  If 'verbose' is an
    integer > 1, further messages from 'scanqtl' are also
    displayed.
}
  \item{tol}{
    Tolerance for convergence for the binary trait model.
}
  \item{maxit.fitqtl}{
    Maximum number of iterations for fitting the binary trait
    model.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A 'lodprofileM2' object.
%% Should write more later
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
'plotlodmatlist','plotlodmatlist2', 'profileLodMatfn2'
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
