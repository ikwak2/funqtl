scantwoF <- function(cross, pheno.cols, usec=c("slod","mlod"), n.perm, ...) {
    
    n = nind(cross)
    usec <- match.arg(usec)
    
    if (!missing(pheno.cols))
    pheno.cols = 1:nphe(cross)
    
    if (!all(pheno.cols %in% 1:nphe(cross)))
    stop("pheno.cols should be in a range of 1 to ", nphe(cross))
    
    if (missing(n.perm))
    n.perm <- 0
    
    if (n.perm > 0 ) {
        cross <- calc.genoprob(cross, step=0)
        temp <- cross
        pheno <- cross$pheno[, pheno.cols]
        
        Slods <- NULL;
        Mlods <- NULL;
        Slod <- NULL;
        Mlod <- NULL;
        SlodsH <- NULL;
        SlodsL <- NULL;
        MlodsH <- NULL;
        MlodsL <- NULL;
        
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
            
            Slod <- c(Slod, max(summary(out3)$one))
            Mlod <- c(Mlod, max(summary(out4)$one))
            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )
            SlodsH <- c(SlodsH, max(summary(out3)$lod.int) )
            SlodsL <- c(SlodsL, max(summary(out3)$lod.fv1) )
            MlodsH <- c(MlodsH, max(summary(out4)$lod.int) )
            MlodsL <- c(MlodsL, max(summary(out4)$lod.fv1) )
        }
        
        return( cbind(Slod,Mlod,Slods,Mlods,SlodsH, SlodsL, MlodsH, MlodsL) )
        
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
