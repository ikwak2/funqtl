\name{plotlod}
\alias{plotlod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot signed LOD matrix.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  We can do QTL analysis for all the phenotypes for the function valued
  trait. We will have a LOD matrix that (i,j)th element represent a LOD
  score of ith phenotype at marker position j. We added one more
  information. We ploted -LOD score if the effect of QTL is negative. 
}
\usage{
plotlod(out, eff, y, gap, off.end=, ncol, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    %%     ~~Describe \code{cross} here~~
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{out}{
    %%     ~~Describe \code{z2} here~~
    A matrix scaneone output using "hk" option.
    ex) out <- scanone(cross, pheno.col=1:10, method="hk")
  }
  \item{eff}{
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
  %%     ~~Describe \code{\dots} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
