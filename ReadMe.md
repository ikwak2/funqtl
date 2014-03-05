R package, funqtl
-----------------

IL-YOUP KWAK <ikwak2@stat.wisc.edu>, with contributions from Karl W Broman.

QTL mapping for function valued trait

This is an R package with add-on functions for the [R/qtl package](http://www.rqtl.org) to deal
with function valued trait on QTL mapping.


#### Installation

You need the `install_github` function in [devtools]() package. So install
and load devtools:

    if(!require(devtools)) install.packages("devtools")
	    library(devtools)

And then, use `install_github` function to install R/funqtl:

    install_github("ikwak2/funqtl")
