all: doc vignettes

doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'

vignettes: inst/doc/funqtl.Rmd inst/doc/funqtl.R inst/doc/funqtl.html

inst/doc/funqtl.Rmd: vignettes/funqtl.Rmd
	cp $< $@

vignettes/funqtl.R: vignettes/funqtl.Rmd
	cd vignettes;R -e 'library(knitr);purl("funqtl.Rmd")'

inst/doc/funqtl.R: vignettes/funqtl.R
	cp $< $@

vignettes/funqtl.html: vignettes/funqtl.Rmd
	cd vignettes;R -e 'library(knitr);knit2html("funqtl.Rmd")'

inst/doc/funqtl.html: vignettes/funqtl.html
	cp $< $@
