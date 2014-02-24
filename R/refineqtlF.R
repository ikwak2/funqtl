refineqtlF <-
function (cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL,
    formula, method = c("imp", "hk"), model = c("normal", "binary"),
    verbose = TRUE, maxit = 10, incl.markers = TRUE, keeplodprofile = TRUE,
    tol = 1e-04, maxit.fitqtl = 1000) {

    method <- match.arg(method)
    model <- match.arg(model)
    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)
    if (!is.null(covar) && !is.data.frame(covar)) {
        if (is.matrix(covar) && is.numeric(covar))
            covar <- as.data.frame(covar, stringsAsFactors = TRUE)
        else stop("covar should be a data.frame")
    }
#    if (LikePheVector(pheno.col, nind(cross), nphe(cross))) {
#        cross$pheno <- cbind(pheno.col, cross$pheno)
#        pheno.col <- 1
#    }
    if (!missing(qtl) && (!missing(chr) || !missing(pos) || !missing(qtl.name)))
        warning("qtl argument is provided, and so chr, pos and qtl.name are ignored.")
    if (missing(qtl) && (missing(chr) || missing(pos)))
        stop("Provide either qtl or both chr and pos.")
    if (!missing(qtl)) {
        chr <- qtl$chr
        pos <- qtl$pos
    }
    else {
        if (missing(qtl.name)) {
            if (method == "imp")
                qtl <- makeqtl(cross, chr = chr, pos = pos, what = "draws")
            else qtl <- makeqtl(cross, chr = chr, pos = pos,
                what = "prob")
        }
        else {
            if (method == "imp")
                qtl <- makeqtl(cross, chr = chr, pos = pos, qtl.name = qtl.name,
                  what = "draws")
            else qtl <- makeqtl(cross, chr = chr, pos = pos,
                qtl.name = qtl.name, what = "prob")
        }
    }
    if (method == "imp") {
        if (!("geno" %in% names(qtl))) {
            if ("prob" %in% names(qtl)) {
                warning("The qtl object doesn't contain imputations; using method=\"hk\".")
                method <- "hk"
            }
            else stop("The qtl object needs to be created with makeqtl with what=\"draws\".")
        }
    }
    else {
        if (!("prob" %in% names(qtl))) {
            if ("geno" %in% names(qtl)) {
                warning("The qtl object doesn't contain QTL genotype probabilities; using method=\"imp\".")
                method <- "imp"
            }
            else stop("The qtl object needs to be created with makeqtl with what=\"prob\".")
        }
    }
    if (!all(chr %in% names(cross$geno)))
        stop("Chr ", paste(unique(chr[!(chr %in% cross$geno)]),
            sep = " "), " not found in cross.")
    if (verbose > 1)
        scanqtl.verbose <- TRUE
    else scanqtl.verbose <- FALSE
    cross <- subset(cross, chr = as.character(unique(chr)))
    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        if (method == "imp")
            qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
                what = "draws")
        else qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }
    if (method == "imp" && dim(qtl$geno)[3] != dim(cross$geno[[1]]$draws)[3]) {
        warning("No. imputations in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what = "draws")
    }
    map <- attr(qtl, "map")
    if (is.null(map))
        stop("Input qtl object should contain the genetic map.")
    mind <- min(sapply(map, function(a) {
        if (is.matrix(a)) a <- a[1, ]
        min(diff(a))
    }))/2
    if (mind <= 0)
        mind <- 1e-06
#    if (length(pheno.col) > 1) {
#        pheno.col <- pheno.col[1]
#        warning("refineqtl can take just one phenotype; only the first will be used")
#    }

#    if (is.character(pheno.col)) {
#        num <- find.pheno(cross, pheno.col)
#        if (is.na(num))
#        stop("Couldn't identify phenotype \"", pheno.col,   "\"")
#        pheno.col <- num
#    }

    if (missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    #
    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    #
    pheno <- as.data.frame(cross$pheno[, pheno.cols], stringsAsFactors = TRUE)

    #
    if (!is.null(covar) && nrow(covar) != nrow(pheno))
        stop("nrow(covar) != no. individuals in cross.")
    if (!is.null(covar))
        phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors = TRUE)
    hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if (all(hasmissing))
        stop("All individuals are missing phenotypes or covariates.")
    if (any(hasmissing)) {
        origcross <- cross
        origqtl <- qtl
        cross <- subset(cross, ind = !hasmissing)
#
        pheno <- pheno[!hasmissing,]

        if (!is.null(covar))
            covar <- covar[!hasmissing, , drop = FALSE]
        if (method == "imp")
            qtl$geno <- qtl$geno[!hasmissing, , , drop = FALSE]
        else qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,
            , drop = FALSE])
        qtl$n.ind <- sum(!hasmissing)
    }
    if (missing(formula)) {
        formula <- paste("y ~", paste(qtl$altname, collapse = "+"))
        if (!is.null(covar))
            formula <- paste(formula, "+", paste(colnames(covar),
                collapse = "+"))
        formula <- as.formula(formula)
    }
    if (!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if (all(is.na(m)))
            covar <- NULL
        else covar <- covar[, !is.na(m), drop = FALSE]
    }

#########

    formula <- qtl:::checkformula(formula, qtl$altname, colnames(covar))


# identify which QTL are in the model formula
    tovary <- sort(qtl:::parseformula(formula, qtl$altname, colnames(covar))$idx.qtl)
    if(length(tovary) != qtl$n.qtl)
        reducedqtl <- dropfromqtl(qtl, index=(1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl

    if (any(1:length(tovary) != tovary)) {
        tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
        terms <- strsplit(tempform, " *\\+ *")[[1]]
        for (j in seq(along = terms)) {
            if (length(grep(":", terms[j])) > 0) {
                temp <- strsplit(terms[j], " *: *")[[1]]
                for (k in seq(along = temp)) {
                  g <- grep("^[Qq][0-9]+$", temp[k])
                  if (length(g) > 0) {
                    num <- as.numeric(substr(temp[k], 2, nchar(temp[k])))
                    temp[k] <- paste("Q", which(tovary == num),
                      sep = "")
                  }
                }
                terms[j] <- paste(temp, collapse = ":")
            }
            else {
                g <- grep("^[Qq][0-9]+$", terms[j])
                if (length(g) > 0) {
                  num <- as.numeric(substr(terms[j], 2, nchar(terms[j])))
                  terms[j] <- paste("Q", which(tovary == num),
                    sep = "")
                }
            }
        }
        formula <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    }
    curpos <- pos[tovary]
    chrnam <- chr[tovary]
    if (verbose)
        cat("pos:", curpos, "\n")
    converged <- FALSE
    oldo <- NULL
    lc <- length(chrnam)
    lastout <- vector("list", length(curpos))
    names(lastout) <- qtl$name[tovary]
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    for (i in 1:maxit) {
        if (keeplodprofile) {
#
            basefit <- NULL
            basefitlod <- NULL

            for(phv in pheno.cols) {
                basefit[[phv]] <- qtl:::fitqtlengine(pheno = pheno[,pheno.cols[phv]], qtl = reducedqtl,
                                   covar = covar, formula = formula, method = method,
                                   model = model, dropone = TRUE, get.ests = FALSE,
                                   run.checks = FALSE, cross.attr = cross.attr,
                                   sexpgm = sexpgm, tol = tol, maxit = maxit.fitqtl)
                basefitlod <- c( basefitlod, basefit[[phv]]$result.full[1,4] )
            }

        }

#
        else {
            basefit <- NULL
            basefitlod <- NULL

            for(phv in pheno.cols) {

                basefit[[phv]] <- qtl:::fitqtlengine(pheno = pheno[,pheno.cols[phv]], qtl = reducedqtl,
                          covar = covar, formula = formula, method = method,
                          model = model, dropone = FALSE, get.ests = FALSE,
                          run.checks = FALSE, cross.attr = cross.attr, sexpgm = sexpgm,
                          tol = tol, maxit = maxit.fitqtl)
                basefitlod <- c( basefitlod, basefit[[phv]]$result.full[1,4] )
            }
        }

#
        if (i == 1) {
            origlod <- curlod <- thisitlod <- mean(basefitlod)
            origpos <- curpos
        }
        if (verbose)
            cat("Iteration", i, "\n")
        o <- sample(lc)
        if (!is.null(oldo))
            while (o[1] != oldo[lc]) o <- sample(lc)
        oldo <- o
        newpos <- curpos
        for (j in o) { # j=4
            otherchr <- chrnam[-j]
            otherpos <- newpos[-j]
            thispos <- as.list(newpos)
            if (any(otherchr == chrnam[j])) {
                linkedpos <- otherpos[otherchr == chr[j]]
                if (any(linkedpos < newpos[j]))
                  low <- max(linkedpos[linkedpos < newpos[j]])
                else low <- -Inf
                if (any(linkedpos > newpos[j]))
                  high <- min(linkedpos[linkedpos > newpos[j]])
                else high <- Inf
                thispos[[j]] <- c(low, high)
            }
            else thispos[[j]] <- c(-Inf, Inf)


##
            out <- scanqtlfn(cross = cross, pheno.cols = pheno.cols,
                       chr = chrnam, pos = thispos, covar = covar, formula = formula,
                       method = method, model = model, incl.markers = incl.markers,
                       verbose = scanqtl.verbose, tol = tol, maxit = maxit.fitqtl)

            lastout[[j]] <- out
            newpos[j] <- as.numeric(strsplit(names(out)[out ==
                max(out)], "@")[[1]][2])
            if (verbose) {
                cat(" Q", j, " pos: ", curpos[j], " -> ", newpos[j],
                  "\n", sep = "")
                #
                cat("    LOD increase: ", round(max(out) - curlod,
                  3), "\n")
            }
            curlod <- max(out)
        }
        if (verbose) {
            cat("all pos:", curpos, "->", newpos, "\n")
            cat("LOD increase at this iteration: ", round(curlod -
                thisitlod, 3), "\n")
        }
        thisitlod <- curlod
        if (max(abs(curpos - newpos)) < mind) {
            converged <- TRUE
            break
        }
        curpos <- newpos
        reducedqtl <- replaceqtl(cross, reducedqtl, seq(length(curpos)),
            reducedqtl$chr, curpos, reducedqtl$name)
    }
    if (verbose) {
        cat("overall pos:", origpos, "->", newpos, "\n")
        cat("LOD increase overall: ", round(curlod - origlod,
            3), "\n")
    }
    if (!converged)
        warning("Didn't converge.")
    g <- grep("^.+@[0-9\\.]+$", qtl$name)
    if (length(g) == length(qtl$name))
        thenames <- NULL
    else thenames <- qtl$name
    if (any(hasmissing)) {
        qtl <- origqtl
        cross <- origcross
    }
    for (j in seq(along = tovary)) qtl <- qtl:::replaceqtl(cross, qtl,
        tovary[j], chrnam[j], newpos[j])
    if (!is.null(thenames))
        qtl$name <- thenames

#####
    if (keeplodprofile) {
#
        dropresult <- basefit[[1]]$result.drop
        if (is.null(dropresult)) {
            if (length(lastout) == 1) {
#
                dropresult <- rbind(c(NA, NA, basefit[[1]]$result.full[1, 4]))
                rownames(dropresult) <- names(lastout)
            }
            else stop("There's a problem: need dropresult, but didn't obtain one.")
        }



        rn <- rownames(dropresult)
        qn <- names(lastout)
        for (i in seq(along = lastout)) {
#
            if (length(lastout) == 1) {
                drprest <- NULL
                for( ii in 1:length(pheno.cols) ) {
                    drprest <- c(drprest, basefit[[ii]]$result.full[1,4] )
                }

                drprest <- mean(drprest)
            } else {
                drprest <- NULL
                for( ii in 1:length(pheno.cols) ) {
                    drprest <- c(drprest, basefit[[ii]]$result.drop[rn == qn[i], 3] )
                }
                drprest <- mean(drprest)
            }
#
            lastout[[i]] <- lastout[[i]] - (max(lastout[[i]]) - drprest)

            pos <- as.numeric(matrix(unlist(strsplit(names(lastout[[i]]),
                "@")), byrow = TRUE, ncol = 2)[, 2])
            chr <- rep(qtl$chr[tovary][i], length(pos))
            lastout[[i]] <- data.frame(chr = chr, pos = pos,
                lod = as.numeric(lastout[[i]]), stringsAsFactors = TRUE)
        }
        names(lastout) <- qtl$name[tovary]
        for (i in seq(along = lastout)) {
            class(lastout[[i]]) <- c("scanone", "data.frame")
            thechr <- qtl$chr[i]
            if (method == "imp")
                detailedmap <- attr(cross$geno[[thechr]]$draws,
                  "map")
            else detailedmap <- attr(cross$geno[[thechr]]$prob,
                "map")
            if (is.matrix(detailedmap))
                detailedmap <- detailedmap[1, ]
            r <- range(lastout[[i]][, 2]) + c(-1e-05, 1e-05)
            rn <- names(detailedmap)[detailedmap >= r[1] & detailedmap <=
                r[2]]
            o <- grep("^loc-*[0-9]+", rn)
            if (length(o) > 0)
                rn[o] <- paste("c", thechr, ".", rn[o], sep = "")
            if (length(rn) == nrow(lastout[[i]]))
                rownames(lastout[[i]]) <- rn
        }
        attr(qtl, "lodprofile") <- lastout
    }
    if ("pLOD" %in% names(attributes(qtl)) && curlod > origlod)
        attr(qtl, "pLOD") <- attr(qtl, "pLOD") + curlod - origlod
    qtl
}
