\name{plotlodmatlist}
\alias{plotlodmatlist}

\title{
  Plot profile LOD maxrix.
}
\description{
  Plot profile LOD maxrix of given multiple QTL model.
}
\usage{
plotlodmatlist(lodmatlist, ylab = "QTL position", xlab = "Time", mval = 0, ...)
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
'plotlodmatlist2', 'profileLodMatfn', 'profileLodMatfn2'
}
\examples{
data(simspal)
simspal <- calc.genoprob(simspal, step=0)
qtlslod <- makeqtl(simspal, chr = c(1, 4),
               pos = c(36.6, 27.8), what = "prob")
lodmat1.c <- profileLodMatfn(simspal, qtl =  qtlslod, pheno.cols =1:241,
                             formula = y~Q1 + Q2 , method = "hk",
                                                          verbose = F)
plotlodmatlist(lodmat1.c, main="The Profile LOD image of data")



# difference between profileLodMatfn1 and 2
# version 1 gives separate two images if there is 2 qtl on same chromosome
# version 2 gives one image for 2 qtl on same chromosome.

par(mfrow=c(1,2))
qtlslod <- makeqtl(simspal, chr = c(1, 1, 4),
               pos = c(36.6, 61, 27.8), what = "prob")
lodmat1 <- profileLodMatfn(simspal, qtl =  qtlslod, pheno.cols =1:241,
                             formula = y~Q1 + Q2 + Q3 , method = "hk")
lodmat2 <- profileLodMatfn2(simspal, qtl =  qtlslod, pheno.cols =1:241,
                             formula = y~Q1 + Q2 + Q3 , method = "hk")

plotlodmatlist(lodmat1, main="The Profile LOD image of data")
plotlodmatlist2(lodmat2, main="The Profile LOD image of data")

}
\keyword{ hplot }