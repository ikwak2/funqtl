\name{plotlod}
\alias{plotlod}
\title{
  Plot signed LOD matrix.
}
\description{
  Plot a signed version of a LOD matrix (taking -LOD when the QTL effect
  is negative), where (i,j) element is for phenotype i at marker j.
}
\usage{
plotlod(output, effects, y, gap=25, off.end=0.5, ncol=251, ...)
}
\arguments{
  \item{output}{
    A matrix scaneone output using "hk" option.
    ex) out <- scanone(cross, pheno.col=1:10, method="hk")
  }
  \item{effects}{
    The sign information whether the QTL having AA do negative effect.
    Get this buy using geteffects() function.
  }
  \item{y}{
    Positions of phenotypes in image.
  }
  
  \item{gap}{
    The gap among chromosomes.
  }

  \item{off.end}{
    The margin between end of chromosome image and a line segment.
  }
  
  \item{ncol}{
    The number of color between red and blue
  }
  
\item{\dots}{
  More graphical components of 'image.plot'.
}
}
\value{
  A graph of signed LOD matrix.
}
\author{
  Karl W Broman, Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>

}

\seealso{
'geteffects'

}
\examples{

data(simspal)
simspal <- calc.genoprob(simspal, step=1)

out <- scanone(simspal, pheno.col=1:241 , method="hk")
eff <- geteffects(simspal, pheno.cols=1:241)

plotlod(out, eff, gap=15)



}
\keyword{ hplot }
