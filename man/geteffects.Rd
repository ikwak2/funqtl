\name{geteffects}
\alias{geteffects}

\title{
  get coefficients of qtl. (get sign information) 
}
\description{
  required function to run plotlod which plot signed lod image. The
  plotlod function need a sign information for each LOD score. This
  geteffects function get coefficients at each gene position. So we can
  get sign information from this function.

}
\usage{
geteffects(cross, pheno.cols)
}
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{pheno.cols}{
    phenotype columns to be used
  }
  
}
\value{
  A matrix of coefficients. (i,j)th item is a coefficient of jth
  position as a qtl of ith
  observation.
}
\author{
  Karl W Broman, Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>

}

\seealso{
'plotlod'

}
\examples{

data(simspal)
simspal <- calc.genoprob(simspal, step=1)

out <- scanone(simspal, pheno.col=1:241 , method="hk")
eff <- geteffects(simspal, pheno.cols=1:241)

plotlod(out, eff, gap=15)



}
\keyword{ models }
