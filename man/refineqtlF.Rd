\name{refineqtlF}
\alias{refineqtlF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Refine the positions of QTL for function valued trait
}
\description{
  Extended version of 'refineqtl' function in 'qtl' package for
  function valued trait.
  Iteratively scan the positions for QTL in the context of a
  multiple QTL model, to try to identify the positions that maximize
  slod criteria, for a fixed QTL model.

}
\usage{
refineqtlF(cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL, formula, method = c("imp", "hk"), model = c("normal", "binary"), verbose = TRUE, maxit = 10, incl.markers = TRUE, keeplodprofile = TRUE, tol = 1e-04, maxit.fitqtl = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
  }
  \item{qtl}{
    A QTL object, as produced by 'makeqtl', containing the
    positions of the QTL.  Provide either 'qtl' or the pair 'chr'
    and 'pos'.
  }
  \item{chr}{
    Vector indicating the chromosome for each QTL; if 'qtl' is
    provided, this should not be.
  }
  \item{pos}{
    Vector indicating the positions for each QTL; if 'qtl' is
    provided, this should not be.
  }
  \item{qtl.name}{
    Optional user-specified name for each QTL.  If 'qtl' is
    provided, this should not be.
  }
  \item{covar}{
    A matrix or data.frame of covariates.  These must be strictly
    numeric.
  }
  \item{formula}{
    An object of class 'formula' indicating the model to be
    fitted.  (It can also be the character string representation
    of a formula.)  QTLs are indicated as 'Q1', 'Q2', etc.
    Covariates are indicated by their names in 'covar'.
  }
  \item{method}{
    Indicates whether to use multiple imputation or Haley-Knott
    regression.
  }
  \item{model}{
    The phenotype model: the usual model or a model for binary
    traits
  }
  \item{verbose}{
    If TRUE, give feedback about progress.  If 'verbose' is an
    integer > 1, further messages from 'scanqtl' are also
    displayed.
  }
  \item{maxit}{
    Maximum number of iterations.

  }
  \item{incl.markers}{
    If FALSE, do calculations only at points on an evenly
    spaced grid.
  }
  \item{keeplodprofile}{
    If TRUE, keep the LOD profiles from the last iteration
    as attributes to the output.
  }
  \item{tol}{
    Tolerance for convergence for the binary trait model.
  }
  \item{maxit.fitqtl}{
    Maximum number of iterations for fitting the binary trait
    model.
  }
}
\details{
  This is an extended version of 'refineqtl' of 'qtl' package. For a
  multiple qtl model, this refines each qtl position to move for a
  better position that fits slod criteria given other qtl positions
  fixed. Do this process iteratively to find a refined version of
  multiple QTL.

}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...



  An object of class 'qtl', with QTL placed in their new positions.


  %%               %%  need to check this !!!  %%
  If 'keeplodprofile=TRUE', LOD profiles from the last pass through
  the refinement algorithm are retained as an attribute,
  '"lodprofile"', to the object.  These may be plotted with
  'plotLodProfile'.
 %%

}
\references{
  %% ~put references to the literature/web site here ~


     Zeng, Z.-B., Kao, C.-H., and Basten, C. J. (1999) Estimating the
     genetic architecture of quantitative traits.  _Genet. Res._ *74*,
     279-289.

     Haley, C. S. and Knott, S. A. (1992) A simple regression method
     for mapping quantitative trait loci in line crosses using flanking
     markers.  _Heredity_ *69*, 315-324.

     Sen, S. and Churchill, G. A. (2001) A statistical framework for
     quantitative trait mapping.  _Genetics_ *159*, 371-387.


}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
     'fitqtl', 'makeqtl', 'scanqtl', 'scanqtlF', 'addtoqtl', 'dropfromqtl',
     'replaceqtl', 'plotLodProfile'
}
\examples{

data(exd)

qtl1.c <- makeqtl(exd, chr = c(2), pos = c(30), what = "prob")
thisqtl1.c <- refineqtlF(exd, pheno.cols = 1:10,  qtl = qtl1.c,
                      formula = y~Q1, method = "hk")


}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
