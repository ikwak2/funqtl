scantwoF <- function(cross, pheno.cols, usec=c("slod","mlod"), n.perm, ...) {

    n = nind(cross)
    usec <- match.arg(usec)

    if (missing(pheno.cols))
    pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
    stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    if (missing(n.perm))
    n.perm <- 0

    if (n.perm > 0 ) {
        temp <- cross
        pheno <- cross$pheno[, pheno.cols]

        Slod <- NULL;
        Mlod <- NULL;
        SlodsH <- NULL;
        SlodsL <- NULL;
        MlodsH <- NULL;
        MlodsL <- NULL;
        Slods <- NULL;
        Mlods <- NULL;

        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]

            out2 <- scantwo(temp, pheno.col = pheno.cols,  ...)
            out1 <- scanone(temp, pheno.col = pheno.cols,  ...)

            # out3 for slod
            out3 <- out2
            out3$lod <- apply(out2$lod, 1:2, mean)

            # out4 for mlod
            out4 <- out2
            out4$lod <- apply(out2$lod, 1:2, max)

            SLOD <- rowMeans(out1[,-(1:2)])
            MLOD <- apply(out1[,-(1:2)], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )
            SlodsH <- c(SlodsH, max(summary(out3)$lod.int) )
            SlodsL <- c(SlodsL, max(summary(out3)$lod.fv1) )
            MlodsH <- c(MlodsH, max(summary(out4)$lod.int) )
            MlodsL <- c(MlodsL, max(summary(out4)$lod.fv1) )
        }


#        out <- list(one_slod = as.matrix(Slods), one_mlod = Mlods, fullvadd_slod = SlodsH, fullvadd_mlod = MlodsH, fv1_slod = SlodsL, fv1_mlod = MlodsL)
        out <- list(one = cbind(Slods, Mlods), fullvadd = cbind(SlodsH, MlodsH), fv1 = cbind(SlodsL, MlodsL) )
        class(out) <- c("scantwoperm", "list")
        out

    } else {
        out <- scantwo(cross, pheno.col = pheno.cols, ...)

        if(usec=="slod") {
            out$lod <- apply(out$lod, 1:2, mean)
        }
        if(usec=="mlod") {
            out$lod <- apply(out$lod, 1:2, max)
        }
        out
    }
}
