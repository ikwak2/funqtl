\name{plot2lod}
\alias{plot2lod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot LOD matrix.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  We can do QTL analysis for all the phenotypes for the function valued
  trait. We will have a LOD matrix that (i,j)th element represent a LOD
  score of ith phenotype at marker position j.
}
\usage{
plot2lod(cross, out1, ylab = "QTL position", xlab = "Time", mval = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
%%     ~~Describe \code{cross} here~~
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{out1}{
%%     ~~Describe \code{z2} here~~
    A scanone output of LOD scores. 
}
  \item{ylab}{
    A lable of y axis. Default is 'QTL position'.
}
  \item{xlab}{
    A lable of x axis. Default is 'Time'.
}
\item{mval}{
  %%     ~~Describe \code{mval} here~~
  The maximum LOD value of legend. The color of legend goes 0 to 'mval'. If this
  value is less than the max lod score, it automatically changed to max
  value.
}
\item{\dots}{
  %%     ~~Describe \code{\dots} here~~
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

plot2lod(simspal,out1, col=heat.colors(100)[100:1],
main="The LOD image of the data")




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
