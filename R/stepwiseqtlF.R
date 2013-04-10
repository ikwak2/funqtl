stepwiseqtlF <- function (cross, chr, pheno.cols, qtl, usec=c("slod","mlod"), formula, max.qtl = 10,
covar = NULL, method = c("imp", "hk"), model = c("normal",
"binary"), incl.markers = TRUE, refine.locations = TRUE,
additive.only = FALSE, penalties, 
keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000)
{
    
    if (!missing(pheno.cols))
    pheno.cols = 1:nphe(cross)
    
    #
    if (!all(pheno.cols %in% 1:nphe(cross)))
    stop("pheno.cols should be in a range of 1 to ", nphe(cross))
    
    
    pheno <- cross$pheno
    
    if (!("cross" %in% class(cross)))
    stop("Input should have class \"cross\".")
    if (!missing(chr))
    cross <- subset(cross, chr)
    #    if (LikePheVector(pheno.col, nind(cross), nphe(cross))) {
    #        cross$pheno <- cbind(pheno.col, cross$pheno)
    #        pheno.col <- 1
    #    }
    if (!missing(qtl)) {
        if (!("qtl" %in% class(qtl)))
        stop("The qtl argument must be an object of class \"qtl\".")
        m <- is.na(match(qtl$chr, names(cross$geno)))
        if (any(m)) {
            wh <- qtl$chr[m]
            if (length(wh) > 1)
            stop("Chromosomes ", paste(wh, collapse = ", "),
            " (in QTL object) not in cross object.")
            else stop("Chromosome ", wh, " (in QTL object) not in cross object.")
        }
        if (missing(formula)) {
            if (!is.null(covar))
            formula <- paste("y ~ ", paste(names(covar),
            collapse = "+"), "+")
            else formula <- "y ~ "
            formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr),
            sep = ""), collapse = "+"))
        }
        else {
            temp <- qtl:::checkStepwiseqtlStart(qtl, formula, covar)
            qtl <- temp$qtl
            formula <- temp$formula
        }
        startatnull <- FALSE
    }
    else {
        if (!missing(formula))
        warning("formula ignored if qtl is not provided.")
        startatnull <- TRUE
    }
    if (!startatnull)
    qtl$name <- qtl$altname
    
    method <- match.arg(method)
    model <- match.arg(model)
    usec <- match.arg(usec)
    
    if (method == "imp") {
        if (!("draws" %in% names(cross$geno[[1]]))) {
            if ("prob" %in% names(cross$geno[[1]])) {
                warning("The cross doesn't contain imputations; using method=\"hk\".")
                method <- "hk"
            }
            else stop("You need to first run sim.geno.")
        }
    }
    else {
        if (!("prob" %in% names(cross$geno[[1]]))) {
            if ("draws" %in% names(cross$geno[[1]])) {
                warning("The cross doesn't contain QTL genotype probabilities; using method=\"imp\".")
                method <- "imp"
            }
            else stop("You need to first run calc.genoprob.")
        }
    }
    if (method == "imp")
    qtlmethod <- "draws"
    else qtlmethod <- "prob"
    if (!missing(qtl) && qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if (method == "imp")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
        what = "draws")
        else qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
        what = "prob")
    }
    if (!missing(qtl) && method == "imp" && dim(qtl$geno)[3] !=
    dim(cross$geno[[1]]$draws)[3]) {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what = "draws")
    }
    if (!startatnull) {
        if (method == "imp" && !("geno" %in% names(qtl)))
        stop("The qtl object doesn't contain imputations; re-run makeqtl with what=\"draws\".")
        else if (method == "hk" && !("prob" %in% names(qtl)))
        stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }
    
    if (!is.null(covar))
    phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors = TRUE)
    hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if (all(hasmissing))
    stop("All individuals are missing phenotypes or covariates.")
    if (any(hasmissing)) {
        pheno <- pheno[!hasmissing]
        cross <- subset(cross, ind = !hasmissing)
        if (!is.null(covar))
        covar <- covar[!hasmissing, , drop = FALSE]
        if (!startatnull) {
            if (method == "imp")
            qtl$geno <- qtl$geno[!hasmissing, , , drop = FALSE]
            else {
                for (i in seq(along = qtl$prob)) qtl$prob[[i]] <- qtl$prob[[i]][!hasmissing,
                , drop = FALSE]
            }
            qtl$n.ind <- sum(!hasmissing)
        }
    }
    if (max.qtl < 1)
    stop("Need max.qtl > 0 if we are to scan for qtl")
    if (is.null(covar)) {
        lod0 <- 0
        if (startatnull)
        firstformula <- y ~ Q1
        else firstformula <- formula
    }
    #####  Need modification
    
    else {
        lod0 <- length(pheno)/2 * log10(sum((pheno - mean(pheno))^2)/sum(lm(pheno ~
        as.matrix(covar))$resid^2))
        if (startatnull)
        firstformula <- as.formula(paste("y~", paste(names(covar),
        collapse = "+"), "+", "Q1"))
        else firstformula <- formula
    }
    cross.type <- class(cross)[1]
    if (missing(penalties)) {
        if (cross.type == "f2") {
            penalties <- c(3.52, 4.28, 2.69)
        }
        else if (cross.type == "bc") {
            penalties <- c(2.69, 2.62, 1.19)
        }
        else stop("No default penalties available for cross type ",
        cross.type)
    }
    else if (length(penalties) != 3) {
        if (length(penalties) == 1) {
            if (additive.only)
            penalties <- c(penalties, Inf, Inf)
            else stop("You must include a penalty for interaction terms.")
        }
        else {
            if (length(penalties) == 2)
            penalties <- penalties[c(1, 2, 2)]
            else {
                warning("penalties should have length 3")
                penalties <- penalties[1:3]
            }
        }
    }
    if (verbose > 2)
    verbose.scan <- TRUE
    else verbose.scan <- FALSE
    curbest <- NULL
    curbestplod <- 0
    if (verbose)
    cat(" -Initial scan\n")
    if (startatnull) {
        if (additive.only || max.qtl == 1 ) {
            out <- scanoneF(cross, pheno.cols = pheno.cols, method = method,
            model = "normal", addcovar = covar)
            if( usec == "slod") {
                lod <- max(out[, 3], na.rm = TRUE)
                curplod <- calc.plod(lod, c(1, 0, 0), penalties = penalties)
                wh <- which(!is.na(out[, 3]) & out[, 3] == lod)
            }
            if( usec == "mlod") {
                lod <- max(out[, 4], na.rm = TRUE)
                curplod <- calc.plod(lod, c(1, 0, 0), penalties = penalties)
                wh <- which(!is.na(out[, 4]) & out[, 4] == lod)
            }
            
            if (length(wh) > 1)
            wh <- sample(wh, 1)
            qtl <- makeqtl(cross, as.character(out[wh, 1]), out[wh,
            2], "Q1", what = qtlmethod)
            formula <- firstformula
            n.qtl <- 1
        }
        else {
            out <- scantwoF(cross, pheno.cols = pheno.cols, usec=usec, method = method,
            model = "normal", incl.markers = incl.markers,
            addcovar = covar, verbose = verbose.scan)
            lod <- out$lod
            lod1 <- max(diag(lod), na.rm = TRUE)
            plod1 <- calc.plod(lod1, c(1, 0, 0), penalties = penalties)
            loda <- max(lod[upper.tri(lod)], na.rm = TRUE)
            ploda <- calc.plod(loda, c(2, 0, 0), penalties = penalties)
            lodf <- max(lod[lower.tri(lod)], na.rm = TRUE)
            plodf <- calc.plod(lodf, c(2, 0, 1), penalties = penalties)
            if (plod1 > ploda && plod1 > plodf) {
                wh <- which(!is.na(diag(lod)) & diag(lod) ==
                lod1)
                if (length(wh) > 1)
                wh <- sample(wh, 1)
                m <- out$map[wh, ]
                qtl <- makeqtl(cross, as.character(m[1, 1]),
                m[1, 2], "Q1", what = qtlmethod)
                formula <- firstformula
                n.qtl <- 1
                lod <- lod1
                curplod <- plod1
            }
            else if (ploda > plodf) {
                temp <- max(out, what = "add")
                if (nrow(temp) > 1)
                temp <- temp[sample(1:nrow(temp), 1), ]
                qtl <- makeqtl(cross, c(as.character(temp[1,
                1]), as.character(temp[1, 2])), c(temp[1, 3],
                temp[1, 4]), c("Q1", "Q2"), what = qtlmethod)
                formula <- as.formula(paste(qtl:::deparseQTLformula(firstformula),
                "+Q2", sep = ""))
                curplod <- ploda
                lod <- loda
                n.qtl <- 2
            }
            else {
                temp <- max(out, what = "full")
                if (nrow(temp) > 1)
                temp <- temp[sample(1:nrow(temp), 1), ]
                qtl <- makeqtl(cross, c(as.character(temp[1,
                1]), as.character(temp[1, 2])), c(temp[1, 3],
                temp[1, 4]), c("Q1", "Q2"), what = qtlmethod)
                formula <- as.formula(paste(qtl:::deparseQTLformula(firstformula),
                "+Q2+Q1:Q2", sep = ""))
                curplod <- plodf
                lod <- lodf
                n.qtl <- 2
            }
        }
    }
    else {
        if (verbose)
        cat(" ---Starting at a model with", length(qtl$chr),
        "QTL\n")
        if (refine.locations) {
            if (verbose)
            cat(" ---Refining positions\n")
            rqtl <- refineqtlF(cross, pheno.cols = pheno.cols, qtl = qtl,
            covar = covar, formula = formula, method = method,
            verbose = verbose.scan, incl.markers = incl.markers,
            keeplodprofile = FALSE)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                cat(" ---  Moved a bit\n")
            }
            qtl <- rqtl
        }
        
        
        res.full = NULL;
        #        qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl, sep = "")
        qtl$name <- qtl$altname
        
        for(ii in pheno.cols) {
            res.full <- c(res.full, fitqtl(cross, pheno.col = ii, qtl, covar = covar, formula = formula,
            method = method, model = model, dropone = FALSE,
            get.ests = FALSE, run.checks = FALSE, tol = tol,
            maxit = maxit)$result.full[1, 4] )
        }
        if(usec=="slod") {
            lod <- mean(res.full) - lod0
        }
        if(usec=="mlod") {
            lod <- max(res.full) - lod0
        }
        
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
        penalties = penalties)
        attr(qtl, "pLOD") <- curplod
        n.qtl <- length(qtl$chr)
    }
    attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod
    if (curplod > 0) {
        curbest <- qtl
        curbestplod <- curplod
        if (verbose)
        cat("** new best ** (pLOD increased by ", round(curplod,
        4), ")\n", sep = "")
    }
    if (keeptrace) {
        temp <- list(chr = qtl$chr, pos = qtl$pos)
        attr(temp, "formula") <- qtl:::deparseQTLformula(formula)
        attr(temp, "pLOD") <- curplod
        class(temp) <- c("compactqtl", "list")
        thetrace <- list(`0` = temp)
    }
    if (verbose)
    cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
    qtl:::deparseQTLformula(formula), "\n")
    if (verbose > 1)
    cat("         qtl:", paste(qtl$chr, round(qtl$pos, 1),
    sep = "@"), "\n")
    i <- 0
    while (n.qtl < max.qtl) {
        i <- i + 1
        if (verbose) {
            cat(" -Step", i, "\n")
            cat(" ---Scanning for additive qtl\n")
        }
        out <- addqtlF(cross, pheno.cols = pheno.cols, qtl = qtl,
        covar = covar, formula = formula, method = method,
        incl.markers = incl.markers, verbose = verbose.scan)
        
        if(usec=="slod") {
            curlod <- max(out[, 3], na.rm = TRUE)
            wh <- which(!is.na(out[, 3]) & out[, 3] == curlod)
        }
        if(usec=="mlod") {
            curlod <- max(out[, 4], na.rm = TRUE)
            wh <- which(!is.na(out[, 4]) & out[, 4] == curlod)
        }
        
        if (length(wh) > 1)
        wh <- sample(wh, 1)
        curqtl <- addtoqtl(cross, qtl, as.character(out[wh, 1]),
        out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
        curformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
        "+Q", n.qtl + 1, sep = ""))
        curlod <- curlod + lod
        curplod <- calc.plod(curlod, qtl:::countqtlterms(curformula,
        ignore.covar = TRUE), penalties = penalties)
        if (verbose)
        cat("        plod =", curplod, "\n")
        curnqtl <- n.qtl + 1
        if (!additive.only) {
            for (j in 1:n.qtl) { #j=2
                if (verbose)
                cat(" ---Scanning for QTL interacting with Q",
                j, "\n", sep = "")
                thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                "+Q", n.qtl + 1, "+Q", j, ":Q", n.qtl + 1,
                sep = ""))
                out <- addqtlF(cross, pheno.cols = pheno.cols, qtl = qtl,
                covar = covar, formula = thisformula, method = method,
                incl.markers = incl.markers, verbose = verbose.scan)
                
                
                if(usec=="slod") {
                    thislod <- max(out[, 3], na.rm = TRUE)
                    wh <- which(!is.na(out[, 3]) & out[, 3] == thislod)
                }
                if(usec=="mlod") {
                    thislod <- max(out[, 4], na.rm = TRUE)
                    wh <- which(!is.na(out[, 4]) & out[, 4] == thislod)
                }
                
                if (length(wh) > 1)
                wh <- sample(wh, 1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,
                1]), out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
                thislod <- thislod + lod
                thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                ignore.covar = TRUE), penalties = penalties)
                if (verbose)
                cat("        plod =", thisplod, "\n")
                if (thisplod > curplod) {
                    curformula <- thisformula
                    curplod <- thisplod
                    curlod <- thislod
                    curqtl <- thisqtl
                    curnqtl <- n.qtl + 1
                }
            }
            if (n.qtl > 1) {
                if (verbose)
                cat(" ---Look for additional interactions\n")
                
                
                ## <<
                
                ## <<
                
                temp <- addint(cross, pheno.col = pheno.cols[1], qtl,
                covar = covar,
                formula = formula, method = method,
                qtl.only = TRUE,
                verbose = verbose.scan)
                if(!is.null(temp)) {
                    
                    lodlod <- NULL;
                    for(ii in pheno.cols) {
                        lodlod <- cbind(lodlod, addint(cross, pheno.col = ii, qtl,
                        covar = covar,
                        formula = formula, method = method,
                        qtl.only = TRUE,
                        verbose = verbose.scan)[,3] )
                    }
                    if(usec=="slod") {
                        if(!(is.matrix(lodlod))) {
                            lodlod <- mean(lodlod)
                        } else {
                            lodlod <- apply(lodlod,1,mean)
                        }
                        thislod <- max(lodlod, na.rm=TRUE)
                    }
                    
                    if(usec=="mlod") {
                        if(!(is.matrix(lodlod))) {
                            lodlod <- max(lodlod)
                        } else {
                            lodlod <- apply(lodlod,1,max)
                        }
                        thislod <- max(lodlod, na.rm=TRUE)
                    }
                    
                    
                    
                    wh <- which(!is.na(lodlod) & lodlod == thislod)
                    if (length(wh) > 1)
                    wh <- sample(wh, 1)
                    thisformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                    "+", rownames(temp)[wh]))
                    thislod <- thislod + lod
                    thisplod <- calc.plod(thislod, qtl:::countqtlterms(thisformula,
                    ignore.covar = TRUE), penalties = penalties)
                    if (verbose)
                    cat("        plod =", thisplod, "\n")
                    if (thisplod > curplod) {
                        curformula <- thisformula
                        curplod <- thisplod
                        curlod <- thislod
                        curqtl <- qtl
                        curnqtl <- n.qtl
                    }
                }
            }
            #            if (scan.pairs) {
            #                if (verbose)
            #                    cat(" ---Scan for an additional pair\n")
            #                out <- addpair(cross, pheno.col = pheno.col,
            #                  qtl = qtl, covar = covar, formula = formula,
            #                  method = method, incl.markers = incl.markers,
            #                  verbose = verbose.scan)
            #                thelod <- out$lod
            #                loda <- max(thelod[upper.tri(thelod)], na.rm = TRUE)
            #                ploda <- calc.plod(loda + lod, c(2, 0, 0, 0) +
            #                  countqtlterms(formula, ignore.covar = TRUE),
            #                  penalties = penalties)
            #                lodf <- max(thelod[lower.tri(thelod)], na.rm = TRUE)
            #                plodf <- calc.plod(lodf + lod, c(2, 0, 1, 1) +
            #                  countqtlterms(formula, ignore.covar = TRUE),
            #                  penalties = penalties)
            #                if (verbose) {
            #                  cat("        ploda =", ploda, "\n")
            #                  cat("        plodf =", plodf, "\n")
            #                }
            #                if (ploda > curplod && loda > plodf) {
            #                  temp <- max(out, what = "add")
            #                  if (nrow(temp) > 1)
            #                    temp <- temp[sample(1:nrow(temp), 1), ]
            #                  curqtl <- addtoqtl(cross, qtl, c(as.character(temp[1,
            #                    1]), as.character(temp[1, 2])), c(temp[1,
            #                    3], temp[1, 4]), paste("Q", n.qtl + 1:2,
            #                    sep = ""))
            #                  curformula <- as.formula(paste(deparseQTLformula(formula),
            #                    "+Q", n.qtl + 1, "+Q", n.qtl + 2, sep = ""))
            #                  curplod <- ploda
            #                  lod <- loda + lod
            #                  curnqtl <- n.qtl + 2
            #                }
            #                else if (plodf > curplod) {
            #                  temp <- max(out, what = "full")
            #                  if (nrow(temp) > 1)
            #                    temp <- temp[sample(1:nrow(temp), 1), ]
            #                  curqtl <- addtoqtl(cross, qtl, c(as.character(temp[1,
            #                    1]), as.character(temp[1, 2])), c(temp[1,
            #                    3], temp[1, 4]), paste("Q", n.qtl + 1:2,
            #                    sep = ""))
            #                  curformula <- as.formula(paste(deparseQTLformula(formula),
            #                    "+Q", n.qtl + 1, "+Q", n.qtl + 2, "+Q", n.qtl +
            #                      1, ":Q", n.qtl + 2, sep = ""))
            #                  curplod <- plodf
            #                  lod <- lodf + lod
            #                  curnqtl <- n.qtl + 2
            #                }
            #            }
        }
        qtl <- curqtl
        n.qtl <- curnqtl
        attr(qtl, "formula") <- qtl:::deparseQTLformula(curformula)
        attr(qtl, "pLOD") <- curplod
        formula <- curformula
        lod <- curlod
        if (refine.locations) {
            if (verbose)
            cat(" ---Refining positions\n")
            rqtl <- refineqtlF(cross, pheno.cols = pheno.cols, qtl = qtl,
            covar = covar, formula = formula, method = method,
            verbose = verbose.scan, incl.markers = incl.markers,
            keeplodprofile = FALSE)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                cat(" ---  Moved a bit\n")
                qtl <- rqtl
                
                
                res.full = NULL;
                for(ii in pheno.cols) {
                    res.full <- c(res.full, fitqtl(cross, pheno.col = ii, qtl,
                    covar = covar, formula = formula,
                    method = method, model = model,
                    dropone = FALSE,
                    get.ests = FALSE, run.checks = FALSE,
                    tol = tol,
                    maxit = maxit)$result.full[1, 4] )
                }
                if(usec=="slod") {
                    lod <- mean(res.full) - lod0
                }
                if(usec=="mlod") {
                    lod <- max(res.full) - lod0
                }
                
                
                curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                ignore.covar = TRUE), penalties = penalties)
                attr(qtl, "pLOD") <- curplod
            }
        }
        
        
        if (verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
        "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
        cat("         qtl:", paste(qtl$chr, round(qtl$pos,
        1), sep = "@"), "\n")
        if (curplod > curbestplod) {
            if (verbose)
            cat("** new best ** (pLOD increased by ", round(curplod -
            curbestplod, 4), ")\n", sep = "")
            curbest <- qtl
            curbestplod <- curplod
        }
        if (keeptrace) {
            temp <- list(chr = qtl$chr, pos = qtl$pos)
            attr(temp, "formula") <- qtl:::deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- i
            thetrace <- c(thetrace, temp)
        }
        if (n.qtl >= max.qtl)
        break
    }
    
    if (verbose)
    cat(" -Starting backward deletion\n")
    while (n.qtl > 1) {
        i <- i + 1
        
        ## <<
        
        #        cat(qtl$name)
        #        cat(qtl$altname)
        #        cat("\n")
        qtl$name <- qtl$altname
        out2 <- fitqtl(cross, pheno.col=pheno.cols[1], qtl, covar = covar, formula = formula,
        method = method, model = model, dropone = TRUE, get.ests = FALSE,
        run.checks = FALSE, tol = tol, maxit = maxit)$result.drop
        rn <- rownames(out2)
        wh <- c(grep("^[Qq][0-9]+$", rn), grep("^[Qq][0-9]+:[Qq][0-9]+$", rn))
        
        ## <<
        outout <- NULL;
        for(ii in pheno.cols) {
            outout <- cbind(outout, fitqtl(cross, pheno.col=ii, qtl, covar = covar,
            formula = formula, method = method,
            model = model, dropone = TRUE, get.ests = FALSE,
            run.checks = FALSE, tol = tol,
            maxit = maxit)$result.drop[, 3]
            )
        }
        
        
        
        if(usec=="slod") {
            outout <- apply(outout,1,mean)
        }
        if(usec=="mlod") {
            outout <- apply(outout,1,max)
        }
        out <- outout[wh , drop = FALSE]
        thelod <- out
        minlod <- min(thelod, na.rm = TRUE)
        
        
        wh <- which(!is.na(thelod) & thelod == minlod)
        if (length(wh) > 1)
        wh <- sample(wh, 1)
        lod <- lod - minlod
        todrop <- rn[wh]
        if (verbose)
        cat(" ---Dropping", todrop, "\n")
        if (length(grep(":", todrop)) > 0) {
            theterms <- attr(terms(formula), "factors")
            wh <- colnames(theterms) == todrop
            if (!any(wh))
            stop("Confusion about what interation to drop!")
            theterms <- colnames(theterms)[!wh]
            formula <- as.formula(paste("y~", paste(theterms,
            collapse = "+")))
        }
        else {
            numtodrop <- as.numeric(substr(todrop, 2, nchar(todrop)))
            theterms <- attr(terms(formula), "factors")
            cn <- colnames(theterms)
            g <- c(grep(paste("^[Qq]", numtodrop, "$", sep = ""),
            cn), grep(paste("^[Qq]", numtodrop, ":", sep = ""),
            cn), grep(paste(":[Qq]", numtodrop, "$", sep = ""),
            cn))
            cn <- cn[-g]
            formula <- as.formula(paste("y~", paste(cn, collapse = "+")))
            if (n.qtl > numtodrop) {
                for (j in (numtodrop + 1):n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
                j, j - 1)
            }
            qtl <- dropfromqtl(qtl, index = numtodrop)
            qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl,
            sep = "")
            n.qtl <- n.qtl - 1
        }
        curplod <- calc.plod(lod, qtl:::countqtlterms(formula, ignore.covar = TRUE),
        penalties = penalties)
        if (verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
        "  formula:", qtl:::deparseQTLformula(formula), "\n")
        if (verbose > 1)
        cat("         qtl:", paste(qtl$chr, round(qtl$pos,
        1), sep = ":"), "\n")
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- curplod
        if (refine.locations) {
            if (verbose)
            cat(" ---Refining positions\n")
            if (!is.null(qtl)) {
                rqtl <- refineqtlF(cross, pheno.cols = pheno.cols,
                qtl = qtl, covar = covar, formula = formula,
                method = method, verbose = verbose.scan, incl.markers = incl.markers,
                keeplodprofile = FALSE)
                if (any(rqtl$pos != qtl$pos)) {
                    if (verbose)
                    cat(" ---  Moved a bit\n")
                    qtl <- rqtl
                    
                    # <<
                    
                    
                    
                    res.full = NULL;
                    for(ii in pheno.cols) {
                        res.full <- c(res.full, fitqtl(cross, pheno.col = ii, qtl,
                        covar = covar, formula = formula,
                        method = method, model = model,
                        dropone = FALSE,
                        get.ests = FALSE, run.checks = FALSE,
                        tol = tol,
                        maxit = maxit)$result.full[1, 4] )
                    }
                    if(usec=="slod") {
                        lod <- mean(res.full) - lod0
                    }
                    if(usec=="mlod") {
                        lod <- max(res.full) - lod0
                    }
                    
                    curplod <- calc.plod(lod, qtl:::countqtlterms(formula,
                    ignore.covar = TRUE), penalties = penalties)
                    attr(qtl, "pLOD") <- curplod
                }
            }
        }
        if (curplod > curbestplod) {
            if (verbose)
            cat("** new best ** (pLOD increased by ", round(curplod -
            curbestplod, 4), ")\n", sep = "")
            curbestplod <- curplod
            curbest <- qtl
        }
        if (keeptrace) {
            temp <- list(chr = qtl$chr, pos = qtl$pos)
            attr(temp, "formula") <- qtl:::deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- i
            thetrace <- c(thetrace, temp)
        }
    }
    
    
    
    
    
    
    if (!is.null(curbest)) {
        chr <- curbest$chr
        pos <- curbest$pos
        o <- order(factor(chr, levels = names(cross$geno)), pos)
        qtl <- makeqtl(cross, chr[o], pos[o], what = qtlmethod)
        formula <- as.formula(attr(curbest, "formula"))
        if (length(chr) > 1) {
            n.qtl <- length(chr)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
            i, n.qtl + i)
            for (i in 1:n.qtl) formula <- qtl:::reviseqtlnuminformula(formula,
            n.qtl + o[i], i)
        }
        #        if (keeplodprofile) {
        #    if (verbose)
        #    cat(" ---One last pass through refineqtl\n")
        #    qtl <- refineqtl(cross, pheno.col = pheno.col, qtl = qtl,
        #    covar = covar, formula = formula, method = method,
        #    verbose = verbose.scan, incl.markers = incl.markers,
        #    keeplodprofile = TRUE)
        #}
        attr(qtl, "formula") <- qtl:::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- attr(curbest, "pLOD")
        curbest <- qtl
    }
    else {
        curbest <- numeric(0)
        class(curbest) <- "qtl"
        attr(curbest, "pLOD") <- 0
    }
    if (keeptrace)
    attr(curbest, "trace") <- thetrace
    attr(curbest, "formula") <- qtl:::deparseQTLformula(attr(curbest,
    "formula"), TRUE)
    curbest
}
