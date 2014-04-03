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
  \code{\link{plotlod}}
}
\examples{

data(simspal)
\dontshow{simspal <- subset(simspal,chr=1:4,ind=1:100)}
# Genotype probabilities for H-K
\dontrun{simspal <- calc.genoprob(simspal, step=2.5)
}\dontshow{simspal <- calc.genoprob(simspal)}

out <- scanone(simspal, pheno.col=100:241 , method="hk")

# Calculate signed information
eff <- geteffects(simspal, pheno.cols=100:241)

# plot signed lod image
plotlod(out, eff, gap=15)



}
\keyword{ models }
