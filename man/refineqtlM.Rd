\name{refineqtlM}
\alias{refineqtlM}

\title{
  Refine the positions of QTL for function valued trait
}
\description{
  Extended version of 'refineqtl' function in 'qtl' package for
  function valued trait.
  Iteratively scan the positions for QTL in the context of a
  multiple QTL model, to try to identify the positions that maximize
  "hk" and "f" criteria, for a fixed QTL model.
}
\usage{
refineqtlM(cross, Y, qtl, chr, pos, qtl.name, formula, verbose = TRUE,
  maxit = 10, incl.markers = TRUE, tol = 1e-04, maxit.fitqtl = 1000, method=c("hk","f"), pheno.cols)
}

\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{Y}{ Matrix of phenotypes }
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
  \item{formula}{
    An object of class 'formula' indicating the model to be
    fitted.  (It can also be the character string representation
    of a formula.)  QTLs are indicated as 'Q1', 'Q2', etc.
    Covariates are indicated by their names in 'covar'.
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
  \item{tol}{
    Tolerance for convergence for the binary trait model.
  }
  \item{maxit.fitqtl}{
    Maximum number of iterations for fitting the binary trait
    model.
  }
  \item{method}{
    Indicates whether to use "hk" or "f" criteria
  }
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
  }
}
\details{
  This is an extended version of 'refineqtl' of 'qtl' package. For a
  multiple qtl model, this refines each qtl position to move for a
  better position that fits "hk" or "f" criteria given other qtl positions
  fixed. Do this process iteratively to find a refined version of
  multiple QTL.

}
\value{
  An object of class 'qtl', with QTL placed in their new positions.
}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}


\seealso{
          \code{\link{refineqtlM}}, \code{\link[qtl]{fitqtl}}, \code{\link[qtl]{makeqtl}},
     \code{\link[qtl]{scanqtl}},
     \code{\link[qtl]{addtoqtl}}, \code{\link[qtl]{dropfromqtl}}, \code{\link[qtl]{replaceqtl}}
}
\examples{

# will be added

}
\keyword{ models }
