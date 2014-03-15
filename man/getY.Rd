\name{getY}
\alias{getY}

\title{
  Do dimensional reduction using pca.
}
\description{
  Do dimensional reduction using pca.
}
\usage{
getY(cross, n.max=5, criteria=.9, nn=0)
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
  'scanoneM' 
}
\examples{

     data(exd)
     
     exd <- calc.genoprob(exd, step=2)
     Y <- getY(exd, criteria=.9)
     out1 <- scanoneM(exd, Y, method = "hk")
     out2 <- scanoneM(exd, Y, method = "f")
     out3 <- scanoneM(exd, Y, method = "sl")
     out4 <- scanoneM(exd, Y, method = "ml")
 

     # Summarize results
     summary(out)
     
     # Plot the results : red for slod, blue for mlod
     plot(out, out[,c(1,2,4)], col = c("red","blue"), ylab = "lod")
     

     # Permutation tests
     perm1 <- scanone(exd, method="hk", n.perm=100)
     # Threshold from the permutation test
     summary(perm1, alpha=c(0.05, 0.10))
     
     # Results above the 0.05 threshold
     summary(out, perms=perm1, alpha=0.05)
     
     
}
\keyword{ utilities }
