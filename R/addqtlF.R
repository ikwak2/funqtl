addqtlF <-
function(cross, pheno.cols, ...) {

    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    LODS <- NULL;
    for(i in pheno.cols ) {
        out <- addqtl(cross, pheno.col = i, ...)
        LODS <- cbind(LODS, out$lod)
    }

    MXy <- max(LODS)
    Slods <- apply(LODS, 1, mean)
    Mlods <- apply(LODS, 1, max)
    out[,3] <- Slods
    out[,4] <- Mlods
    names(out)[3:4] <- c("slod","mlod")

    out[,1:4]
}
