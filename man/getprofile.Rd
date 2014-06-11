\name{getprofile}
\alias{getprofile}
\title{Make a profile LOD score matrix for a given multiple QTL model.}
\usage{
getprofile(cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL, formula,
  method = c("hk", "imp"), verbose = TRUE, tpy = c("comb", "sep"))
}
\arguments{
  \item{cross}{An object of class \code{"cross"}. See
  \code{\link[qtl]{read.cross}} for detail.}

  \item{pheno.cols}{Columns in the phenotype matrix to be
  used as the phenotype.}

  \item{qtl}{A QTL object, as produced by
  \code{\link[qtl]{makeqtl}}, containing the positions of
  the QTL.  Provide either \code{qtl} or the pair
  \code{chr} and \code{pos}.}

  \item{chr}{Vector indicating the chromosome for each QTL;
  if \code{qtl} is provided, this should not be.}

  \item{pos}{Vector indicating the positions for each QTL;
  if \code{qtl} is provided, this should not be.}

  \item{qtl.name}{Optional user-specified name for each
  QTL.  If \code{qtl} is provided, this should not be.}

  \item{covar}{A matrix or data.frame of covariates.  These
  must be strictly numeric.}

  \item{formula}{An object of class 'formula' indicating
  the model to be fitted.  (It can also be the character
  string representation of a formula.) QTLs are indicated
  as 'Q1', 'Q2', etc.  Covariates are indicated by their
  names in 'covar'.}

  \item{method}{Indicates whether to use multiple
  imputation or Haley-Knott regression.}

  \item{verbose}{If TRUE, give feedback about progress.  If
  'verbose' is an integer > 1, further messages from
  'scanqtl' are also displayed.}

  \item{tpy}{type of output. If there are more QTL's in one
  chromosome. It plot them separately if \code{tpy =
  "sep"}, On the other hand, it combine then in one
  chromosome taking maximum values of them if \code{tpy =
  "comb"}.}
}
\value{
A \code{"lodprofileM"} or \code{"lodprofileM2"} object.
}
\description{
Calculate profile LOD scores for each chromosome.
}
\examples{
data(simspal)
simspal <- calc.genoprob(simspal, step=0)
\dontshow{simspal <- subset(simspal,chr=c(1,3,4), ind=1:50)}

phe <- 1:nphe(simspal)
\dontshow{phe <- seq(1, nphe(simspal), by=60)}

# difference between tpy = "sep" and "comb"

par(mfrow=c(1,2))
qtlslod <- makeqtl(simspal, chr = c(1, 1, 4),
               pos = c(36.6, 61, 27.8), what = "prob")
lodmat1 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =phe,
                    formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "sep")
plotprofile(lodmat1, main="tpy=\\"sep\\"")

lodmat2 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =phe,
                    formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "comb")
plotprofile(lodmat2, main="tpy=\\"comb\\"")
}
\author{
Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}
\seealso{
\code{\link{plotprofile}}
}
\keyword{models}

