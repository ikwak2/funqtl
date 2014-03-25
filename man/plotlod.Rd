\name{plotlod}
\alias{plotlod}
\title{
  Plot signed or plain LOD image.
}
\description{
  Plot signed or plain LOD image. (taking -LOD when the QTL effect
  is negative), where (i,j) element is for phenotype i at marker j.
}
\usage{
plotlod(output, effects, y, gap=25, ncol=251, rev = FALSE, ...)
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

  \item{off.end}{
    ..
  }

  \item{gap}{
    The gap among chromosomes.
  }

  \item{ncol}{
    The number of color between red and blue
  }

  \item{rev}{
    transepose x and y axis
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
  \code{\link{geteffects}}

}
\examples{

data(simspal)
simspal <- calc.genoprob(simspal, step=1)

out <- scanone(simspal, pheno.col=1:241 , method="hk")
eff <- geteffects(simspal, pheno.cols=1:241)

# set time increment
nam <- names(out)[-(1:2)]
y <- as.numeric(substr(nam, 2, nchar(nam)))/60

# plot with effects
plotlod(out, eff, y, gap=15)

# without effects and rev = T
plotlod(out, y, gap=15, rev = T)

# plot without y
plotlod(out, eff, gap=25)


}
\keyword{ hplot }
