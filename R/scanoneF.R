scanoneF <-
function(cross, pheno.cols, n.perm, ...) {

    n = nind(cross)

    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    if (missing(n.perm))
        n.perm <- 0

    if (n.perm > 0 ) {

        cross <- calc.genoprob(cross, step=0)
        temp <- cross
        pheno <- cross$pheno[,pheno.cols]

        Slods <- NULL;
        Mlods <- NULL;
        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]

            out <- scanone(temp, pheno.col = pheno.cols, ...)
            SLOD <- rowMeans(out[,-(1:2)])
            MLOD <- apply(out[,-(1:2)], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )

        }
        return( rbind(Slods,Mlods) )

    } else {

        out <- scanone(cross, pheno.col = pheno.cols, ...)
        SLOD <- rowMeans(out[,-(1:2)])
        MLOD <- apply(out[,-(1:2)], 1, max)

        out[,3] <- SLOD
        out[,4] <- MLOD
        names(out)[3:4] <- c("slod","mlod")

        out[,1:4]
    }
}
