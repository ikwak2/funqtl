extmat <- function(matlist, schr, pheno.cols, ...) {
    for (ch in schr) {
        matlist[[as.character(ch)]]$nqtl = 0
        
        ou <- NULL
        k <- 1
        for (i in pheno.cols) {
            outo <- addqtl(chr = ch, pheno.col = i, ...)
            if(k == 1) {
                nn <- nrow(outo)
                ou <- cbind(rep(ch,nn), outo$pos, outo$lod)
                rownames(ou) <- rownames(outo)
                colnames(ou) <- colnames(outo)
                k = 2
            } else {
                ou <- cbind(ou, outo$lod)
            }
        }
        matlist[[as.character(ch)]]$out <- ou
    }
    #    matlist
    chrs <- as.numeric(lapply(matlist, function(t) t$out[1,1]))
    oo <- order(chrs)
    
    matt <- NULL
    for( i in oo) {
        matt[[chrs[i]]] <- matlist[[i]]
    }
    matt
}