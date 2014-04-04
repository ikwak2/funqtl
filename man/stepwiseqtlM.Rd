\name{stepwiseqtlM}
\alias{stepwiseqtlM}
\title{
  Stepwise selection for multiple QTL in function valued trait data
}
\description{
Extension of 'stepwiseqtl' function of 'qtl' package.
Performs forward/backward selection to identify a multiple QTL
model for function valued trait data, with model choice made via a
penalized LOD score, with separate penalties on main effects and interactions.

}
\usage{
stepwiseqtlM(cross, chr, Y,  qtl,
formula, max.qtl = 10, incl.markers=TRUE, refine.locations = TRUE,
    penalties, additive.only = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000,
 method=c("hk", "f"),
 pheno.cols)

}
\arguments{
  \item{cross}{
    An object of class 'cross'. See 'read.cross' for details.
  }
  \item{chr}{
    Optional vector indicating the chromosomes to consider in
    search for QTL.  This should be a vector of character strings
    referring to chromosomes by name; numeric values are
    converted to strings.  Refer to chromosomes with a preceding
    '-' to have all chromosomes but those considered.  A logical
    (TRUE/FALSE) vector may also be used.
  }
  \item{Y}{ Matrix of phenotypes }
  \item{pheno.cols}{
    Columns in the phenotype matrix to be used as the phenotype.
  }
  \item{qtl}{
    Optional QTL object (of class '"qtl"', as created by
    'makeqtl') to use as a starting point.
  }
  \item{formula}{
    Optional formula to define the QTL model to be used as a
    starting point.
  }
  \item{max.qtl}{
    Maximum number of QTL to which forward selection should
    proceed.
  }

  \item{method}{
    Indicates whether to use multiple imputation or Haley-Knott
    regression.
  }

  \item{incl.markers}{
    If FALSE, do calculations only at points on an evenly
    spaced grid.
  }
  \item{refine.locations}{
    If TRUE, use 'refineqtl' to refine the QTL locations
    after each step of forward and backward selection.
  }
  \item{additive.only}{
    If TRUE, allow only additive QTL models; if FALSE,
    consider also pairwise interactions among QTL.
  }
  \item{penalties}{
    Vector of three values indicating the penalty on main
    effects and heavy and light penalties on interactions.  See
    the Details below. If missing, default values are used that
    are based on simulations of backcrosses and intercrosses with
    genomes modeled after that of the mouse.
  }
  \item{keeptrace}{
    If TRUE, keep information on the sequence of models visited
    through the course of forward and backward selection as an
    attribute to the output.
  }
  \item{verbose}{
    If TRUE, give feedback about progress.  If 'verbose' is an
    integer > 1, even more information is printed.
  }
  \item{tol}{
    Tolerance for convergence for the binary trait model.
  }
  \item{maxit}{
    Maximum number of iterations for fitting the binary trait
    model.
  }

}
\value{

  The output is a representation of the best model, as measured by
  the penalized LOD score (see Details), among all models visited.
  This is QTL object (of class '"qtl"', as produced by 'makeqtl'),
  with attributes '"formula"', indicating the model formula, and
  '"pLOD"' indicating the penalized LOD score.

%%  If 'keeplodprofile=TRUE', LOD profiles from the last pass through
%%  the refinement algorithm are retained as an attribute,
%%  '"lodprofile"', to the object.  These may be plotted with
%%  'plotLodProfile'.

  If 'keeptrace=TRUE', the output will contain an attribute
  '"trace"' containing information on the best model at each step of
  forward and backward elimination.  This is a list of objects of
  class '"compactqtl"', which is similar to a QTL object (as
  produced by 'makeqtl') but containing just a vector of chromosome
     IDs and positions for the QTL.  Each will also have attributes
     '"formula"' (containing the model formula) and '"pLOD"'
     (containing the penalized LOD score.
     If 'n.perm' is missing, the function returns a data.frame whose
     first two columns contain the chromosome IDs and cM positions.
     Subsequent third and fourth columns contain the SLOD and MLOD scores.

     If 'n.perm' is specified, the function returns the results of a
     permutation test and the output returns the matrix of two
     columns. The first column for SLOD and the second column for MLOD score.

}
\author{
  Il-Youp Kwak, <email: ikwak2@stat.wisc.edu>
}


\seealso{
  \code{\link[qtl]{stepwiseqtl}}, \code{\link{stepwiseqtlF}}
}
\examples{

# I'm working on it

}
