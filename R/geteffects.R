geteffects <- function(cross,pheno.cols) {
    if(missing(pheno.cols))
        pheno.cols=1:nphe(cross)

#    out <- scanone(cross, pheno.col = pheno.cols, method="hk")
    phe <- as.matrix(cross$pheno)
    eff <- NULL
    for(i in 1:nchr(cross)) {
        pr <- cross$geno[[i]]$prob[,,2]
        eff <- rbind(eff, t(apply(pr, 2, function(a,b) lm(b~a)$coef[2,], phe)))
    }
    eff
}
