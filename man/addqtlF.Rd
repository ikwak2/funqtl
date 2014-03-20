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
An object of class 'scanone', as produced by the 'scanone'
function.  LOD scores are relative to the base model (with any
terms that include the new QTL omitted).

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


\seealso{
     \code[qtl]{\link{addqtl}}, \code[qtl]{\link{scanone}}, \code[funqtl]{\link{scanoneF}}, \code[qtl]{\link{makeqtl}}
}
\examples{
data(simspal)
qtl1 <- makeqtl(simspal, chr = c(4),
               pos = c(27.8), what = "prob")

added <- addqtlF(simspal, qtl =  qtl1, pheno.cols =1:241,
                             formula = y~Q1 + Q2, method = "hk")

summary(added)



}
\keyword{ models }
