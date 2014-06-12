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
  method <- match.arg(method)

  # drop individuals with missing data
  pheno <- cross$pheno[,pheno.cols, drop=FALSE]
  if(!is.null(covar))
    phcovar <- cbind(pheno, covar)
  else phcovar <- as.data.frame(pheno, stringsAsFactors = TRUE)
  hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
  if (all(hasmissing))
    stop("All individuals are missing phenotypes or covariates.")
  if (any(hasmissing)) {
    pheno <- pheno[!hasmissing,,drop=FALSE]
    cross <- subset(cross, ind = !hasmissing)
    if (!is.null(covar))
      covar <- covar[!hasmissing, , drop = FALSE]
    if (method == "imp")
      qtl$geno <- qtl$geno[!hasmissing, , , drop = FALSE]
    else {
      for (i in seq(along = qtl$prob)) qtl$prob[[i]] <- qtl$prob[[i]][!hasmissing,
                      , drop = FALSE]
    }
    qtl$n.ind <- sum(!hasmissing)
  }

  if(missing(lod0)) {
    if(is.null(covar)) lod0 <- rep(0, length(pheno.cols))
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
    lod[i] <- fitqtl(cross=cross, pheno.col=pheno.cols[i], qtl=qtl, covar=covar, formula=formula,
                     method=method, model="normal", dropone=FALSE, get.ests=FALSE,
                     run.checks=FALSE)$result.full[1,4]
  
  lod - lod0                   
}
