\name{plotlodmatlist2}
\alias{plotlodmatlist2}

\title{
  Plot profile LOD maxrix. (Version 2)
}
\description{
  Plot profile LOD maxrix of given multiple QTL model.
}
\usage{
plotlodmatlist2(lodmatlist, ylab = "QTL position", xlab = "Time", mval = 0, ...)
}
\arguments{
  \item{lodmatlist}{
    An 'lodprofileM2' object that produced by 'profileLodMatfn2' function.
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
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

\seealso{
'plotlodmatlist', 'profileLodMatfn', 'profileLodMatfn2'
}
\examples{
}
\keyword{ hplot }
