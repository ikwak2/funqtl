\name{plot2lod}
\alias{plot2lod}

\title{
  Plot LOD matrix.
}
\description{
  Plot a matrix of LOD scores, for which the (i,j)th element is
  for phenotype i at marker j.
}
\usage{
plot2lod(cross, out1, ylab = "QTL position", xlab = "Time", mval = 0,
         col=rev(heat.colors(100)), ...)
}

\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{out1}{
    A scanone output of LOD scores. 
}
  \item{ylab}{
    A lable of y axis. Default is 'QTL position'.
}
  \item{xlab}{
    A lable of x axis. Default is 'Time'.
}
\item{mval}{
  The maximum LOD value of legend. The color of legend goes 0 to 'mval'. If this
  value is less than the max lod score, it automatically changed to max
  value.
}
\item{col}{
  Vector of colors
}
\item{\dots}{
  More graphical components of 'image.plot'.
}
}
\value{
  A graph of LOD matrix.
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>

}

\seealso{
'plotlodmatlist', 'plotlodmatlist2'

}
\examples{

# call dataset
data(simspal)
simspal <- calc.genoprob(simspal, step=1)
out1 <- scanone(simspal, pheno.col=1:241 , method="hk")

\dontshow{par(mar=c(3.1, 3.1, 1.6, 0.6))}
plot2lod(simspal, out1, main="The LOD image of the data")




}
\keyword{ hplot }
