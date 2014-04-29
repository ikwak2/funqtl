# fitqtlF: fit a multiple-qtl model to each of several phenotypes
#          return vector of LOD scores
#
# @param cross cross object
# @param pheno.cols
# @param qtl qtl object
# @param formula QTL formula
# @param covar covariate matrix (data.frame)
# @param method whether to use haley-knott regression or imputation to fit the model
# @param lod0 null log10 likelihood for each phenotype (optional): a vector of length ncol(pheno)
fitqtlF <-
function(cross, pheno.cols, qtl, formula, covar=NULL, method=c("hk", "imp"), lod0)
{
  if(missing(lod0)) {
    if(is.null(covar)) lod0 <- 0
    else {
      pheno <- cross$pheno[,pheno.cols,drop=FALSE]
      rss0 <- colSums(lm(as.matrix(pheno) ~ as.matrix(covar))$resid^2, na.rm=TRUE)
      rss00 <- colSums(lm(as.matrix(pheno) ~ 1)$resid^2, na.rm=TRUE)
      lod0 <- nrow(pheno)/2 * log10(rss00/rss0)
    }
  }

  stopifnot(length(lod0) == length(pheno.cols))
  
  lod <- rep(NA, length(pheno.cols))
  for(i in seq(along=pheno.cols))
    lod[i] <- fitqtl(cross=cross, pheno.col=i, qtl=qtl, covar=covar, formula=formula,
                     method=method, model="normal", dropone=FALSE, get.ests=FALSE,
                     run.checks=FALSE)$result.full[1,4]
  
  lod - lod0                   
}
