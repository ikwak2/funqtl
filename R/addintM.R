addintM <- function (cross, Y, qtl, formula, qtl.only = FALSE, verbose = TRUE,
    pvalues = TRUE, simple = FALSE, tol = 1e-04, maxit = 1000, method=c("hk","f"), pheno.cols=pheno.cols)
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
    if (qtl$n.ind != nind(cross)) {
        warning("No. individuals in qtl object doesn't match that in the input cross; re-creating qtl object.")
        qtl <- makeqtl(cross, qtl$chr, qtl$pos, qtl$name,
            what = "prob")
    }

    n.qtl <- qtl$n.qtl
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                sep = "+")
        formula <- as.formula(formula)
    }
    formula <- qtl::checkformula(formula, qtl$altname, NULL)
    factors <- attr(terms(formula), "factors")
    if (sum(factors[1, ]) == 0)
        factors <- factors[-1, ]
    fn <- fn.alt <- rownames(factors)
    qan <- qtl$altname
    qn <- qtl$name
    m <- match(fn, qan)
    fn.alt[!is.na(m)] <- qn[m[!is.na(m)]]
    int2test <- int2test.alt <- NULL
    for (i in 1:(nrow(factors) - 1)) {
        for (j in (i + 1):nrow(factors)) {
            temp <- rep(0, nrow(factors))
            temp[c(i, j)] <- 1
            if (!any(apply(factors, 2, function(a, b) all(a ==
                b), temp))) {
                int2test <- c(int2test, paste(fn[i], fn[j], sep = ":"))
                int2test.alt <- c(int2test.alt, paste(fn.alt[i],
                  fn.alt[j], sep = ":"))
            }
        }
    }
    if (qtl.only && length(int2test) > 0) {
        z <- matrix(unlist(strsplit(int2test, ":")), ncol = 2,
            byrow = TRUE)
        wh <- apply(z, 1, function(a) length(grep("^[Qq][0-9]+$",
            a)))
        int2test <- int2test[wh == 2]
        int2test.alt <- int2test.alt[wh == 2]
    }
    n2test <- length(int2test)
    if (n2test == 0) {
        if (verbose)
            cat("No pairwise interactions to add.\n")
        return(NULL)
    }
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)

###
    thefit0 <- fitqtlM(cross = cross, Y = Y, qtl = qtl,
        formula = formula, tol = tol, method=method, pheno.cols=pheno.cols)

    results <- matrix(ncol = 7, nrow = n2test)
    dimnames(results) <- list(int2test.alt, c("df", "Type III SS",
        "LOD", "%var", "F value", "Pvalue(Chi2)", "Pvalue(F)"))
    for (k in seq(along = int2test)) {


        thefit1 <- fitqtlM(cross = cross, Y = Y, qtl = qtl,
            formula = as.formula(paste(qtl::deparseQTLformula(formula),
                int2test[k], sep = "+")), tol = tol, method=method, pheno.cols=pheno.cols)

        results[k, 1] <- thefit1$result.full[1, 1] - thefit0$result.full[1,
            1]
        results[k, 2] <- thefit1$result.full[1, 2] - thefit0$result.full[1,
            2]
        results[k, 3] <- thefit1$result.full[1, 4] - thefit0$result.full[1,
            4]
        results[k, 4] <- 100 * (1 - 10^(-2 * thefit1$result.full[1,
            4]/qtl$n.ind)) - 100 * (1 - 10^(-2 * thefit0$result.full[1,
            4]/qtl$n.ind))
        results[k, 5] <- (results[k, 2]/results[k, 1])/thefit1$result.full[2,
            3]
        results[k, 6] <- pchisq(results[k, 3] * 2 * log(10),
            results[k, 1], lower.tail = FALSE)
        results[k, 7] <- pf(results[k, 5], results[k, 1], thefit1$result.full[3,
            1], lower.tail = FALSE)
    }
    results <- as.data.frame(results, stringsAsFactors = TRUE)
    class(results) <- c("addint", "data.frame")
    attr(results, "method") <- "hk"
    attr(results, "model") <- "normal"
    attr(results, "formula") <- qtl::deparseQTLformula(formula)
    if (simple)
        pvalues <- FALSE
    attr(results, "pvalues") <- pvalues
    attr(results, "simple") <- simple
    results
}




