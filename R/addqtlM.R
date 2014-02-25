addqtlM <- function (cross, Y, chr, qtl, formula,
    incl.markers = TRUE, verbose = FALSE, tol = 1e-04, maxit = 1000, method = c("hk","f", pheno.cols)
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
    if (!("qtl" %in% class(qtl)))
        stop("The qtl argument must be an object of class \"qtl\".")
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if (verbose > 1) {
        verbose <- TRUE
        verbose.scanqtl <- TRUE
    } else verbose.scanqtl <- FALSE

    n.qtl <- qtl$n.qtl
    qtlchr <- qtl$chr
    qtlpos <- qtl$pos

    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }
    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
        attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
        stepwidth.var <- TRUE
        incl.markers <- TRUE
    } else stepwidth.var <- FALSE


    if (missing(chr)) {
        chr <- names(cross$geno)
    } else chr <- qtl:::matchchr(chr, names(cross$geno))


    n.covar <- 0
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                                sep = "+")
        newformula <- as.formula(paste(formula, "+Q", n.qtl +
            1, sep = ""))
        formula <- as.formula(formula)
    }
    else {
        newqtl <- paste("Q", n.qtl + 1, sep = "")
        formula <- qtl:::checkformula(formula, c(qtl$altname, newqtl),
                                NULL)
        theterms <- rownames(attr(terms(formula), "factors"))
        g <- grep(paste("^[Qq]", n.qtl + 1, "$", sep = ""), theterms)
        if (length(g) == 0) {
            newformula <- as.formula(paste(qtl:::deparseQTLformula(formula),
                                           "+ Q", n.qtl + 1, sep = ""))
        }     else {
            newformula <- formula
            theterms <- colnames(attr(terms(formula), "factors"))
            g <- unique(c(grep(paste("^[Qq]", n.qtl + 1, "$",
                sep = ""), theterms), grep(paste("^[Qq]", n.qtl +
                1, " *:", sep = ""), theterms), grep(paste(": +[Qq]",
                n.qtl + 1, " *:", sep = ""), theterms), grep(paste(": +[Qq]",
                n.qtl + 1, "$", sep = ""), theterms)))
            if (length(g) > 0) {
                theterms <- theterms[-g]
                formula <- as.formula(paste("y ~ ", paste(theterms,
                  collapse = " + "), sep = ""))
            }
        }
    }
    thefactors <- rownames(attr(terms(formula), "factors"))
    todrop <- NULL
    for (i in 1:n.qtl) {
        if (length(grep(paste("^[Qq]", i, "$", sep = ""), thefactors)) ==
            0)
            todrop <- c(todrop, i)
    }
    if (length(todrop) > 0) {
        newqtlnum <- n.qtl + 1
        notdropped <- (1:n.qtl)[-todrop]
        newnum <- 1:length(notdropped)
        qtl <- dropfromqtl(qtl, index = todrop)
        qtlchr <- qtlchr[-todrop]
        qtlpos <- qtlpos[-todrop]
        n.qtl <- n.qtl - length(todrop)
        revnewqtlnum <- n.qtl + 1
        formula <- qtl:::reviseqtlnuminformula(formula, notdropped,
            newnum)
        newformula <- qtl:::reviseqtlnuminformula(newformula, c(notdropped,
            newqtlnum), c(newnum, revnewqtlnum))
    }

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    lod0 <- as.numeric(getlodM(cross, Y, qtl = qtl, formula = formula, tol= tol, method=method, pheno.cols=pheno.cols))



#lod0 <- qtl:::fitqtlengine(pheno = pheno[,1], qtl = qtl, covar = NULL,
#        formula = formula, method = "hk", model = "normal", dropone = FALSE,
#        get.ests = FALSE, run.checks = FALSE, cross.attr = cross.attr,
#        sexpgm = sexpgm, tol = tol, maxit = maxit)$result.full[1, 4]



    results <- NULL
    for (i in chr) {
        if (verbose)
            cat("Scanning chr", i, "\n")
        thechr <- c(qtlchr, i)
        thepos <- c(as.list(qtlpos), list(c(-Inf, Inf)))
        sqout <- scanqtlM(cross, Y, chr = thechr,
            pos = thepos, formula = newformula,
            incl.markers = incl.markers,
            verbose = verbose.scanqtl, tol = tol, maxit = maxit, method=method, pheno.cols=pheno.cols)

        if ("map" %in% names(attributes(cross$geno[[i]]$prob)))
            map <- attr(cross$geno[[i]]$prob, "map")
        else {
            stp <- attr(cross$geno[[i]]$prob, "step")
            oe <- attr(cross$geno[[i]]$prob, "off.end")
            if ("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
                stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
            else stpw <- "fixed"
            map <- create.map(cross$geno[[i]]$map, stp, oe,
                              stpw)
        }

        if (is.matrix(map))
            map <- map[1, ]
        step <- attr(cross$geno[[i]]$prob, "step")
        if (!incl.markers && step > 0) {
            eq.sp.pos <- seq(min(map), max(map), by = step)
            wh.eq.pos <- match(eq.sp.pos, map)
            map <- map[wh.eq.pos]
        }
        w <- names(map)
        o <- grep("^loc-*[0-9]+", w)
        if (length(o) > 0)
            w[o] <- paste("c", i, ".", w[o], sep = "")
        z <- data.frame(lod = as.numeric(sqout) - lod0, stringsAsFactors = TRUE)
        z <- cbind(chr = rep(i, length(map)), pos = as.numeric(map),
            z)
        rownames(z) <- w
        results <- rbind(results, z)
    }
    class(results) <- c("scanone", "data.frame")
    attr(results, "method") <- "hk"
    attr(results, "formula") <- qtl:::deparseQTLformula(newformula)
    results
}






