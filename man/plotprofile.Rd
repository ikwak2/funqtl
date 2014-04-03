\name{plotprofile}
\alias{plotprofile}

\title{
  Plot profile LOD maxrix.
}
\description{
  Plot profile LOD maxrix of given multiple QTL model.
}
\usage{
plotprofile(lodmatlist, ylab = "QTL position", xlab = "Time", mval = 0, ...)
}
\arguments{
  \item{lodmatlist}{
An 'lodprofileM' object that produced by 'profileLodMatfn' function.
}
\item{ylab}{
  A lable of y. Default is "QTL position".
}
  \item{xlab}{
    A lable of x. Default is "Time".
}
  \item{mval}{
  The maximum LOD value of legend. The color of legend goes 0 to 'mval'. If this
  value is less than the max lod score, it automatically changed to max
  value.

}
  \item{\dots}{
More graphical components of 'image.plot'.

}
}
\value{
A graph of profile LOD matrix.
}
\author{
  Il-Youp Kwak, <ikwak2@stat.wisc.edu>
}

\seealso{
  \code{\link{getprofile}}
}
\examples{
data(simspal)
\dontshow{simspal <- subset(simspal,chr=c(1,3,4), ind=1:50)}
# Genotype probabilities for H-K
# it should be step = 0 to run this function
simspal <- calc.genoprob(simspal, step=0)

# difference between tpy = "sep" and "comb"
par(mfrow=c(1,2))
qtlslod <- makeqtl(simspal, chr = c(1, 1, 4),
               pos = c(36.6, 61, 27.8), what = "prob")
lodmat1 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =1:241,
                    formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "sep")
lodmat2 <- getprofile(simspal, qtl =  qtlslod, pheno.cols =1:241,
                    formula = y~Q1 + Q2 + Q3 , method = "hk", tpy = "comb")

plotprofile(lodmat1, main="The Profile LOD image of data")
plotprofile(lodmat2, main="The Profile LOD image of data")

}
\keyword{ hplot }