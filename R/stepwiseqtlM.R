stepwiseqtlM <- function (cross, chr, Y, qtl, formula, max.qtl = 10,
    incl.markers = TRUE, refine.locations = TRUE,
    penalties,  additive.only = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000, method=method, pheno.cols=pheno.cols)
{
    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno)
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")
    if (!missing(chr))
        cross <- subset(cross, chr)

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
            temp <- qtl::checkStepwiseqtlStart(qtl, formula, covar)
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

    qtlmethod <- "prob"
    if (!missing(qtl) && qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }

    if (!startatnull) {
        if (method == "hk" && !("prob" %in% names(qtl)))
            stop("The qtl object doesn't contain QTL genotype probabilities; re-run makeqtl with what=\"prob\".")
    }

    if (max.qtl < 1)
        stop("Need max.qtl > 0 if we are to scan for qtl")
    lod0 <- 0
    if (startatnull)
        firstformula <- y ~ Q1
    else firstformula <- formula

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
            stop("You must include a penalty for interaction terms.")
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
        {
            out <- scanoneM(cross, Y)
            lod <- max(out[, 3], na.rm = TRUE)
            curplod <- calc.plod(lod, c(1, 0, 0), penalties = penalties)
            wh <- which(!is.na(out[, 3]) & out[, 3] == lod)
            if (length(wh) > 1)
                wh <- sample(wh, 1)
            qtl <- makeqtl(cross, as.character(out[wh, 1]), out[wh,
                2], "Q1", what = qtlmethod)
            formula <- firstformula
            n.qtl <- 1
        }

    }  else {
        if (verbose)
            cat(" ---Starting at a model with", length(qtl$chr),
                "QTL\n")
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
            }
            qtl <- rqtl
        }
        lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
            tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
        curplod <- calc.plod(lod, qtl::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        attr(qtl, "pLOD") <- curplod
        n.qtl <- length(qtl$chr)
    }
    attr(qtl, "formula") <- qtl::deparseQTLformula(formula)
    attr(qtl, "pLOD") <- curplod
    if (curplod > 0) {
        curbest <- qtl
        curbestplod <- curplod
        if (verbose)
            cat("** new best ** (pLOD increased by ", round(curplod,
                4), ")\n", sep = "")
    }

    if (verbose)
        cat("    no.qtl = ", n.qtl, "  pLOD =", curplod, "  formula:",
            qtl::deparseQTLformula(formula), "\n")
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
        out <- addqtlM(cross, Y = Y, qtl = qtl,
            formula = formula, incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
        curlod <- max(out[, 3], na.rm = TRUE)
        wh <- which(!is.na(out[, 3]) & out[, 3] == curlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        curqtl <- addtoqtl(cross, qtl, as.character(out[wh, 1]),
            out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
        curformula <- as.formula(paste(qtl::deparseQTLformula(formula),
            "+Q", n.qtl + 1, sep = ""))
        curlod <- curlod + lod
        curplod <- calc.plod(curlod, qtl::countqtlterms(curformula,
            ignore.covar = TRUE), penalties = penalties)
        if (verbose)
            cat("        plod =", curplod, "\n")
        curnqtl <- n.qtl + 1
        if (!additive.only) {
            for (j in 1:n.qtl) {
                if (verbose)
                  cat(" ---Scanning for QTL interacting with Q",
                    j, "\n", sep = "")
                thisformula <- as.formula(paste(qtl::deparseQTLformula(formula),
                  "+Q", n.qtl + 1, "+Q", j, ":Q", n.qtl + 1,
                  sep = ""))
                out <- addqtlM(cross, Y = Y, qtl = qtl,
                  formula = thisformula,
                  incl.markers = incl.markers, verbose = verbose.scan, method=method, pheno.cols=pheno.cols)
                thislod <- max(out[, 3], na.rm = TRUE)
                wh <- which(!is.na(out[, 3]) & out[, 3] == thislod)
                if (length(wh) > 1)
                  wh <- sample(wh, 1)
                thisqtl <- addtoqtl(cross, qtl, as.character(out[wh,
                  1]), out[wh, 2], paste("Q", n.qtl + 1, sep = ""))
                thislod <- thislod + lod
                thisplod <- calc.plod(thislod, qtl::countqtlterms(thisformula,
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



                temp <- addintM(cross, Y = Y, qtl,
                  formula = formula, qtl.only = TRUE,
                  verbose = verbose.scan, method=method, pheno.cols=pheno.cols)






                if (!is.null(temp)) {
                  thislod <- max(temp[, 3], na.rm = TRUE)
                  wh <- which(!is.na(temp[, 3]) & temp[, 3] ==
                    thislod)
                  if (length(wh) > 1)
                    wh <- sample(wh, 1)
                  thisformula <- as.formula(paste(qtl::deparseQTLformula(formula),
                    "+", rownames(temp)[wh]))
                  thislod <- thislod + lod
                  thisplod <- calc.plod(thislod, qtl::countqtlterms(thisformula,
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

        }
        qtl <- curqtl
        n.qtl <- curnqtl
        attr(qtl, "formula") <- qtl::deparseQTLformula(curformula)
        attr(qtl, "pLOD") <- curplod
        formula <- curformula
        lod <- curlod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            rqtl <- refineqtlM(cross, Y = Y, qtl = qtl,
                formula = formula,
                verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)
            if (any(rqtl$pos != qtl$pos)) {
                if (verbose)
                  cat(" ---  Moved a bit\n")
                qtl <- rqtl
                lod <- as.numeric(getlodM(cross, Y = Y, qtl, formula = formula,
                  tol = tol, method=method, pheno.cols=pheno.cols)) - lod0
                curplod <- calc.plod(lod, qtl::countqtlterms(formula,
                  ignore.covar = TRUE), penalties = penalties)
                attr(qtl, "pLOD") <- curplod
            }
        }
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl::deparseQTLformula(formula), "\n")
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
        if (n.qtl >= max.qtl)
            break
    }
    if (verbose)
        cat(" -Starting backward deletion\n")
    while (n.qtl > 1) {
        i <- i + 1

###
        out <- fitqtlM(cross, Y = Y, qtl = qtl, formula = formula, method=method, pheno.cols=pheno.cols)$result.drop


        rn <- names(out)
        wh <- c(grep("^[Qq][0-9]+$", rn), grep("^[Qq][0-9]+:[Qq][0-9]+$",
            rn))
        thelod <- out[wh]
        minlod <- min(thelod, na.rm = TRUE)
        wh <- which(!is.na(thelod) & thelod == minlod)
        if (length(wh) > 1)
            wh <- sample(wh, 1)
        lod <- lod - minlod
        todrop <- names(out)[wh]
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
                for (j in (numtodrop + 1):n.qtl) formula <- qtl::reviseqtlnuminformula(formula,
                  j, j - 1)
            }
            qtl <- dropfromqtl(qtl, index = numtodrop)
            qtl$name <- qtl$altname <- paste("Q", 1:qtl$n.qtl,
                sep = "")
            n.qtl <- n.qtl - 1
        }
        curplod <- calc.plod(lod, qtl::countqtlterms(formula, ignore.covar = TRUE),
            penalties = penalties)
        if (verbose)
            cat("    no.qtl = ", n.qtl, "  pLOD =", curplod,
                "  formula:", qtl::deparseQTLformula(formula), "\n")
        if (verbose > 1)
            cat("         qtl:", paste(qtl$chr, round(qtl$pos,
                1), sep = ":"), "\n")
        attr(qtl, "formula") <- qtl::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- curplod
        if (refine.locations) {
            if (verbose)
                cat(" ---Refining positions\n")
            if (!is.null(qtl)) {


                rqtl <- refineqtlM(cross, Y = Y,
                  qtl = qtl, formula = formula,
                  verbose = verbose.scan, incl.markers = incl.markers, method=method, pheno.cols=pheno.cols)

                if (any(rqtl$pos != qtl$pos)) {
                  if (verbose)
                    cat(" ---  Moved a bit\n")
                  qtl <- rqtl
                  lod <- as.numeric( getlodM(cross, Y = Y, qtl = qtl,
                    formula = formula, tol = tol, method=method, pheno.cols=pheno.cols) ) - lod0
                  curplod <- calc.plod(lod, qtl::countqtlterms(formula,
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
    }
    if (!is.null(curbest)) {
        chr <- curbest$chr
        pos <- curbest$pos
        o <- order(factor(chr, levels = names(cross$geno)), pos)
        qtl <- makeqtl(cross, chr[o], pos[o], what = qtlmethod)
        formula <- as.formula(attr(curbest, "formula"))
        if (length(chr) > 1) {
            n.qtl <- length(chr)
            for (i in 1:n.qtl) formula <- qtl::reviseqtlnuminformula(formula,
                i, n.qtl + i)
            for (i in 1:n.qtl) formula <- qtl::reviseqtlnuminformula(formula,
                n.qtl + o[i], i)
        }
        attr(qtl, "formula") <- qtl::deparseQTLformula(formula)
        attr(qtl, "pLOD") <- attr(curbest, "pLOD")
        curbest <- qtl
    }
    else {
        curbest <- numeric(0)
        class(curbest) <- "qtl"
        attr(curbest, "pLOD") <- 0
    }
    attr(curbest, "formula") <- qtl::deparseQTLformula(attr(curbest,
        "formula"), TRUE)
    curbest
}


