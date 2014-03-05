R package, funqtl
-----------------

IL-YOUP KWAK <ikwak2@stat.wisc.edu>, with contributions from Karl W Broman.

QTL mapping for function valued trait

This is an R package with add-on functions for the [R/qtl package](http://www.rqtl.org) to deal
with function valued trait on QTL mapping. We can do QTL analysis with
function valued trait with this add-on on "qtl" package(made by Karl W
Broman).  

#### Installation

You first need to install [R/qtl](http://www.rqtl.org) 

    if(!require(qtl)) install.packages("qtl")

You also need the `install_github` function in
[Hadley Wickham](http://had.co.nz/)'s [devtools]() package. So install
and load devtools:

    if(!require(devtools)) install.packages("devtools")
	    library(devtools)

Finally, use `install_github` function to install R/funqtl:

    install_github("ikwak2/funqtl")
