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
plotthis <- geteffects(cross)
plotlod(plotthis$out, plotthis$eff, y, gap, off.end=, ncol, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    %%     ~~Describe \code{cross} here~~
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{plotthis$out}{
    %%     ~~Describe \code{z2} here~~
    A matrix scaneone output using "hk" option.
    ex) out <- scanone(cross, pheno.col=1:10, method="hk")
  }
  \item{plotthis$eff}{
    The sign information whether the QTL having AA do negative effect. 
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
'plotlodmatlist', 'plotlodmatlist2'

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cross, z2, ylab = "QTL position", xlab = "Time", mval = 0,
    ...)
{
    nch <- nchr(cross)
    nchs <- lapply(cross$geno, function(t) ncol(t$data))
    end <- NULL
    start <- NULL
    num = 0
    z3 <- NULL
    for (i in 1:nch) {
        z3 <- cbind(z3, cbind(z2[, ((num + 1):(num + nchs[[i]]))],
            rep(mval, nrow(z2)), rep(mval, nrow(z2))))
        num = num + nchs[[i]]
        end <- c(end, num + (i - 1) * 2)
        start <- c(start, num - nchs[[i]] + 1 + (i - 1) * 2)
    }
    midpt <- (start + end)/2
    x <- 1:(dim(z3)[1])
    y <- 1:(dim(z3)[2])
    par(mar = c(5.1, 6.1, 2.1, 2.1))
    image.plot(x, y, z3, yaxt = "n", xlab = xlab, ylab = ylab,
        ...)
    last <- end[nch]
    start <- c(start, last + 2)
    for (i in 1:nch) rect(0, end[i], nrow(z2) + 1, start[i +
        1], col = "white")
    rect(0, end[i], nrow(z2) + 1, start[i + 1] + 1, col = "white")
    u <- par("usr")
    width <- 0.01 * diff(u[1:2])
    text(u[1] - 5 * width, midpt, 1:nch, xpd = TRUE)
    for (i in seq(along = start)) {
        if (i\%\%2)
            rect(u[1] - width, start[i], u[1] - width * 2, end[i],
                col = "gray30", xpd = TRUE)
        else rect(u[1] - width * 2, start[i], u[1] - width *
            3, end[i], col = "gray30", xpd = TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
