\name{plotlodmatlist2}
\alias{plotlodmatlist2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Plot profile LOD maxrix. (Version 2)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Plot profile LOD maxrix of given multiple QTL model.
}
\usage{
plotlodmatlist2(lodmatlist, ylab = "QTL position", xlab = "Time", mval = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lodmatlist}{
    An 'lodprofileM2' object that produced by 'profileLodMatfn2' function.
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
A lable of y. Default is "QTL position".
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
A lable of x. Default is "Time".
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A graph of profile LOD matrix.
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
'plotlodmatlist', 'profileLodMatfn', 'profileLodMatfn2'
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
