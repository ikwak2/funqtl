R package, funqtl
-----------------

IL-YOUP KWAK <ikwak2@stat.wisc.edu>, with contributions from Karl W Broman.

QTL mapping for function valued trait

This is an R package with add-on functions for the [R/qtl package](http://www.rqtl.org) to deal
with function valued trait on QTL mapping.


#### Installation

You need a latest version of [qtl](http://www.rqtl.org) package

    if(!require(devtools)) install.packages("devtools")
    install_github("qtl", "kbroman", "devel")

or for Mac, download
[This](http://www.biostat.wisc.edu/~kbroman/tmp/qtl_1.31-6.tgz) and

    R CMD INSTALL qtl_1.31-6.tgz

And then, use install R/funqtl:

    install_github("ikwak2/funqtl")
