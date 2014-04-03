\name{calcpca}
\alias{calcpca}

\title{
  Do dimensional reduction using pca.
}
\description{
  Do dimensional reduction using pca.
}
\usage{
calcpca(cross, n.max=5, criteria=.9, nn=0)
}

\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
}
  \item{n.max}{
    The number of maximum reduced dimension.
}
  \item{criteria}{
    how much of variance explained.
}
  \item{nn}{
    The number of exact reduced dimension
}

}
\value{
     It gives a matrix that each column have principal components.
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}


\seealso{
  \code{\link{scanoneM}}
}
\examples{
data(simspal)

\dontshow{simspal <- subset(simspal,chr=1:4,ind=1:100)}
# Genotype probabilities for H-K
\dontrun{simspal <- calc.genoprob(simspal, step=2.5)
}\dontshow{simspal <- calc.genoprob(simspal)}

# dimensional reduction of Y
Y <- calcpca(simspal, criteria=.999)

# do multitrait mapping
out.hk <- scanoneM(simspal, Y=Y, method="hk")
out.f  <- scanoneM(simspal, Y=Y, method="f")
out.sl <- scanoneM(simspal, Y=Y, method="sl")
out.ml <- scanoneM(simspal, Y=Y, method="ml")

# Summarize results
summary(out.hk)
summary(out.f)
summary(out.sl)
summary(out.ml)


# Plot the results
par(mfrow=c(3,1))
plot(out.hk)
plot(out.f)
plot(out.sl, out.ml)

}
\keyword{ utilities }
