doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'

build:
	R -e 'library(devtools);build("../funqtl")'

check:
	cd ..; R CMD check funqtl;

vig:
	cd vignettes; R -e 'library(knitr);knit2html("funqtl.Rmd")'
