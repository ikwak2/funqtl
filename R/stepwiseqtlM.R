#' Stepwise selection for multiple QTL in function valued trait data
#'
#'
#' Extension of the R/qtl function \code{\link[qtl]{stepwiseqtl}}. Performs
#' forward/backward selection to identify a multiple QTL model for function
#' valued trait data, with model choice made via a penalized LOD score, with
#' separate penalties on main effects and interactions.
#'
#'
#' @param cross An object of class \code{"cross"}. See \code{\link[qtl]{read.cross}} for details.
#' @param chr Optional vector indicating the chromosomes to consider in search
#' for QTL.  This should be a vector of character strings referring to
#' chromosomes by name; numeric values are converted to strings.  Refer to
#' chromosomes with a preceding \code{"-"} to have all chromosomes but those
#' considered.  A logical (TRUE/FALSE) vector may also be used.
#' @param pheno.cols Columns in the phenotype matrix to be used as the
#' phenotype.
#' @param Y Demension reduced data set. getY(cross) get reduced data set using
#' PCA.
#' @param method which criteria to use: \code{"hk"}, \code{"f"}, \code{"sl"}, or \code{"ml"}.
#' @param qtl Optional QTL object (of class \code{"qtl"}, as created by \code{\link[qtl]{makeqtl}})
#' to use as a starting point.
#' @param formula Optional formula to define the QTL model to be used as a
#' starting point.
#' @param max.qtl Maximum number of QTL to which forward selection should
#' proceed.
#' @param incl.markers If FALSE, do calculations only at points on an evenly
#' spaced grid.
#' @param refine.locations If TRUE, use 'refineqtl' to refine the QTL locations
#' after each step of forward and backward selection.
#' @param additive.only If TRUE, allow only additive QTL models; if FALSE,
#' consider also pairwise interactions among QTL.
#' @param penalties Vector of three values indicating the penalty on main
#' effects and heavy and light penalties on interactions.  See the Details
#' below. If missing, default values are used that are based on simulations of
#' backcrosses and intercrosses with genomes modeled after that of the mouse.
#' @param keeptrace If TRUE, keep information on the sequence of models visited
#' through the course of forward and backward selection as an attribute to the
#' output.
#' @param verbose If TRUE, give feedback about progress.  If 'verbose' is an
#' integer > 1, even more information is printed.
#' @param tol Tolerance for convergence for the binary trait model.
#' @param maxit Maximum number of iterations for fitting the binary trait
#' model.
#' @export
#' @return
#'
#' The output is a representation of the best model, as measured by the
#' penalized LOD score (see Details), among all models visited.  This is QTL
#' object (of class \code{"qtl"}, as produced by \code{\link[qtl]{makeqtl}}), with attributes
#' \code{"formula"}, indicating the model formula, and \code{"pLOD"} indicating the
#' penalized LOD score.
#'
#' If \code{keeptrace=TRUE}, the output will contain an attribute \code{"trace"}
#' containing information on the best model at each step of forward and
#' backward elimination.  This is a list of objects of class \code{"compactqtl"},
#' which is similar to a QTL object (as produced by \code{\link[qtl]{makeqtl}}) but containing
#' just a vector of chromosome IDs and positions for the QTL.  Each will also
#' have attributes \code{"formula"} (containing the model formula) and \code{"pLOD"}
#' (containing the penalized LOD score.  
#'
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{refineqtlF}}, \code{\link{addqtlF}}
#' @references Manichaikul, A., Moon, J. Y., Sen, S, Yandell, B. S. and Broman,
#' K. W. (2009) A model selection approach for the identification of
#' quantitative trait loci in experimental crosses, allowing epistasis.
#' _Genetics_, *181*, 1077-1086.
#'
#' Broman, K. W. and Speed, T. P. (2002) A model selection approach for the
#' identification of quantitative trait loci in experimental crosses (with
#' discussion). _J Roy Stat Soc B_ *64*, 641-656, 731-775.
#'
#' Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
#' quantitative trait loci in line crosses using flanking markers.  _Heredity_
#' *69*, 315-324.
#'
#' Sen, S. and Churchill, G. A. (2001) A statistical framework for quantitative
#' trait mapping.  _Genetics_ *159*, 371-387.
#'
#' Zeng, Z.-B., Kao, C.-H. and Basten, C. J. (1999) Estimating the genetic
#' architecture of quantitative traits.  _Genetical Research_, *74*, 279-289.
#' @examples
#' cat("An example needs to be added.\n")

stepwiseqtlM <- function (cross, chr, Y, qtl, formula, max.qtl = 10,
    incl.markers = TRUE, refine.locations = TRUE,
    penalties,  additive.only = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000, method=c("hk", "f", "sl", "ml"), pheno.cols)
{
    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }

    method <- match.arg(method)

    if(missing(Y)) {
        p <- nphe(cross)
        Y <- as.matrix(cross$pheno[,pheno.cols])
    } else {
        if(is.vector(Y)) { p = 1} else {p = ncol(Y)}
    }

    if (!("cross" %in% class(cross)))
        stop("Input should have class \"cross\".")

    if (method == "sl" || method =="ml" ) {
        temp <- cross
        temp$pheno[, p] <- Y
        if (method == "sl") {
            mtd = "slod"
        } else {
            mtd = "mlod"
        }
        out <- stepwiseqtlF(cross = temp, chr = chr, qtl = qtl,
                            pheno.cols = 1:p, usec = mtd,
                            method = "hk", formula = formula, max.qtl = max.qtl,
                            covar = covar, incl.markers= incl.markers,
                            refine.locations = refine.locations, penalties = penalties,
                            keeptrace = keeptrace, verbose = verbose, tol = tol,
                            maxit = maxit)
        return(out)
        
    }
           
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
            formula <- "y ~ "
            formula <- paste(formula, paste(paste("Q", 1:length(qtl$chr),
                sep = ""), collapse = "+"))
        }
        else {
            temp <- qtl::checkStepwiseqtlStart(qtl, formula)
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
        stop("No default penalties available for cross type ", cross.type)
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


