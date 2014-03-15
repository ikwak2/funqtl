refineqtlM <- function (cross, Y, qtl, chr, pos, qtl.name, formula,
    verbose = TRUE, maxit = 10, incl.markers = TRUE,
    tol = 1e-04, maxit.fitqtl = 1000, method=c("hk","f"), pheno.cols=pheno.cols)
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
        stop("The cross argument must be an object of class \"cross\".")
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)


    if (!missing(qtl) && (!missing(chr) || !missing(pos) || !missing(qtl.name)))
        warning("qtl argument is provided, and so chr, pos and qtl.name are ignored.")
    if (missing(qtl) && (missing(chr) || missing(pos)))
        stop("Provide either qtl or both chr and pos.")
    if (!missing(qtl)) {
        chr <- qtl$chr
        pos <- qtl$pos
    }  else {
        if (missing(qtl.name)) {
            qtl <- makeqtl(cross, chr = chr, pos = pos, what = "prob")
        }
        else {
            qtl <- makeqtl(cross, chr = chr, pos = pos,
                           qtl.name = qtl.name, what = "prob")
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
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name, what = "prob")
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

    if (missing(formula)) {
        formula <- paste("y ~", paste(qtl$altname, collapse = "+"))
        formula <- as.formula(formula)
    }
    formula <- qtl:::checkformula(formula, qtl$altname,NULL)

    tovary <- sort(qtl:::parseformula(formula, qtl$altname, NULL)$idx.qtl)
    if (length(tovary) != qtl$n.qtl)
        reducedqtl <- qtl::dropfromqtl(qtl, index = (1:qtl$n.qtl)[-tovary])
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

        basefit <- getlodM(cross=cross, Y=Y, formula=formula, qtl = qtl, method=method, pheno.cols=pheno.cols)

        if (i == 1) {
            origlod <- curlod <- thisitlod <- as.numeric(basefit)
            origpos <- curpos
        }
        if (verbose)
            cat("Iteration", i, "\n")
        o <- sample(lc)
        if (!is.null(oldo))
            while (o[1] != oldo[lc]) o <- sample(lc)
        oldo <- o
        newpos <- curpos
        for (j in o) {
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
            out <- scanqtlM(cross = cross, Y,
                            chr = chrnam, pos = thispos, formula = formula,
                            incl.markers = incl.markers,
                            verbose = scanqtl.verbose, tol = tol, maxit = maxit.fitqtl)
            lastout[[j]] <- out
            newpos[j] <- as.numeric(strsplit(names(out)[out ==
                                                        max(out)], "@")[[1]][2])
            if (verbose) {
                cat(" Q", j, " pos: ", curpos[j], " -> ", newpos[j],
                    "\n", sep = "")
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

    for (j in seq(along = tovary)) qtl <- replaceqtl(cross, qtl,
                  tovary[j], chrnam[j], newpos[j])
    if (!is.null(thenames))
        qtl$name <- thenames

    if ("pLOD" %in% names(attributes(qtl)) && curlod > origlod)
        attr(qtl, "pLOD") <- attr(qtl, "pLOD") + curlod - origlod
    qtl
}





