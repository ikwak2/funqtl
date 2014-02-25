scanqtlM <- function (cross, Y, chr, pos, formula, incl.markers = FALSE,
                      verbose = TRUE, tol = 1e-04, maxit = 1000, method=c("hk","f"), pheno.cols)
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


    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")

    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)

    if (!("prob" %in% names(cross$geno[[1]]))) {
        stop("You need to first run calc.genoprob.")
    }

    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
        attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
        stepwidth.var <- TRUE
        incl.markers <- TRUE
    } else stepwidth.var <- FALSE

    type <- class(cross)[1]
    chrtype <- sapply(cross$geno, class)
    if (length(chr) != length(pos))
        stop("Input chr and pos must have the same length")
    ichr <- match(chr, names(cross$geno))
    if (any(is.na(ichr)))
        stop("There's no chromosome number ", chr[is.na(ichr)],
            " in input cross object")
    n.qtl <- length(chr)
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                sep = "+")
        formula <- as.formula(formula)
    }
    else {
        formula.str <- qtl:::deparseQTLformula(formula)
        for (i in 1:n.qtl) {
            qtl.term <- paste("Q", i, sep = "")
            if (length(grep(qtl.term, formula.str, ignore.case = TRUE)) ==
                0)
                formula.str <- paste(formula.str, qtl.term, sep = "+")
        }
        formula <- as.formula(formula.str)
    }
    formula <- qtl:::checkformula(formula, paste("Q", 1:length(chr),
        sep = ""), NULL)


    sexpgm <- getsex(cross)
    idx.varied <- NULL
    indices <- pos

##

    for (i in 1:length(pos)) {
        l <- length(pos[[i]])
        if (l >= 2) {
            if (l > 2) {
                msg <- "There are more than two elements in "
                msg <- paste(msg, i, "th input pos.")
                msg <- paste(msg, "The first two are taken as starting and ending position.")
                warning(msg)
            }
            idx.varied <- c(idx.varied, i)

            if ("map" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                map <- attr(cross$geno[[ichr[i]]]$prob, "map")
            else {
                stp <- attr(cross$geno[[ichr[i]]]$prob, "step")
                oe <- attr(cross$geno[[ichr[i]]]$prob, "off.end")
                if ("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$prob)))
                    stpw <- attr(cross$geno[[ichr[i]]]$prob,
                                 "stepwidth")
                else stpw <- "fixed"
                map <- create.map(cross$geno[[ichr[i]]]$map,
                                  stp, oe, stpw)

            }
            if (is.matrix(map))
                map <- map[1, ]
            indices[[i]] <- seq(along = map)
            step <- attr(cross$geno[[ichr[i]]]$prob, "step")
            if (!incl.markers && step > 0) {
                eq.sp.pos <- seq(min(map), max(map), by = step)
                wh.eq.pos <- match(eq.sp.pos, map)
                map <- map[wh.eq.pos]
                indices[[i]] <- indices[[i]][wh.eq.pos]
            }
            start <- pos[[i]][1]
            end <- pos[[i]][2]
            tmp <- which((map - start) <= 0)
            if (length(tmp) != 0)
                start <- map[max(tmp)]
            tmp <- which((end - map) <= 0)
            if (length(tmp) != 0)
                end <- map[min(tmp)]
            pos[[i]] <- as.vector(map[(map >= start) & (map <=
                end)])
            indices[[i]] <- indices[[i]][(map >= start) & (map <=
                end)]
        }
    }
    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)
    n.idx.varied <- length(idx.varied)
    n.loop <- 1
    if (n.idx.varied != 0) {
        idx.pos <- rep(0, n.idx.varied)
        l.varied <- NULL
        for (i in 1:n.idx.varied) {
            l.varied[i] <- length(pos[[idx.varied[i]]])
            n.loop <- n.loop * l.varied[i]
        }
        result <- array(rep(0, n.loop), rev(l.varied))
    }  else {
        qtl <- makeqtl(cross, chr = chr, pos = unlist(pos), what = "prob")

        result <- getlodM(cross, Y, formula=formula, qtl = qtl, tol = tol, method = method, pheno.cols=pheno.cols)



        result <- as.numeric(result)
        names(result) <- "LOD"
        class(result) <- "scanqtl"
        attr(result, "method") <- "hk"
        attr(result, "formula") <- qtl:::deparseQTLformula(formula)
        return(result)
    }
    if (verbose) {
        cat(" ", n.loop, "models to fit\n")
        n.prnt <- floor(n.loop/20)
        if (n.prnt < 1)
            n.prnt <- 1
    }
    current.pos <- NULL
    for (i in 1:n.loop) {
        remain <- i
        if (n.idx.varied > 1) {
            for (j in 1:(n.idx.varied - 1)) {
                ns <- 1
                for (k in (j + 1):n.idx.varied) ns <- ns * length(pos[[idx.varied[k]]])
                idx.pos[j] <- floor(remain/ns) + 1
                remain <- remain - (idx.pos[j] - 1) * ns
                if (remain == 0) {
                  idx.pos[j] <- idx.pos[j] - 1
                  remain <- remain + ns
                }
            }
        }
        idx.pos[n.idx.varied] <- remain
        pos.tmp <- NULL
        for (j in 1:length(pos)) {
            if (j %in% idx.varied) {
                idx.tmp <- which(j == idx.varied)
                pos.tmp <- c(pos.tmp, pos[[j]][idx.pos[idx.tmp]])
            }
            else pos.tmp <- c(pos.tmp, pos[[j]])
        }
        if (is.null(current.pos)) {
            qtl.obj <- makeqtl(cross, chr, pos.tmp, what = "prob")
            current.pos <- pos.tmp
        }
        else {
            thew <- rep(NA, length(pos.tmp))
            for (kk in seq(along = pos.tmp)) {
                if (pos.tmp[kk] != current.pos[kk]) {
                  u <- abs(pos.tmp[kk] - pos[[kk]])
                  w <- indices[[kk]][u == min(u)]
                  if (length(w) > 1) {
                    warning("Confused about QTL positions.  You should probably run jittermap to ensure that no two markers conincide.")
                    w <- sample(w, 1)
                }
                  qtl.obj$prob[[kk]] <- cross$geno[[ichr[kk]]]$prob[, w, ]
                  thew[kk] <- w
                  if (chrtype[ichr[kk]] == "X" && (type == "bc" || type == "f2")) {

                      temp <- qtl.obj$prob[[kk]]
                      temp <- array(temp, dim = c(nrow(temp), 1, ncol(temp)))
                      dimnames(temp) <- list(NULL, "loc", 1:ncol(qtl.obj$prob[[kk]]))
                      qtl.obj$prob[[kk]] <- qtl:::reviseXdata(type,
                      "full", sexpgm, prob = temp, cross.attr = attributes(cross))[, 1, ]
                  }
                  current.pos[kk] <- pos.tmp[kk]
              }
            }
        }

        qtl.obj <- makeqtl(cross, chr, current.pos, what = "prob")
        fit <- getlodM(cross, Y, formula=formula, qtl = qtl.obj, tol=tol, method=method, pheno.cols=pheno.cols)

        if (verbose && ((i - 1)%%n.prnt) == 0)
            cat("    ", i, "/", n.loop, "\n")
        result[i] <- as.numeric(fit)
    }
    dnames <- list(NULL)
    for (i in 1:n.idx.varied) {
        i.chr <- chr[idx.varied[n.idx.varied - i + 1]]
        i.pos <- pos[[idx.varied[n.idx.varied - i + 1]]]
        dnames[[i]] <- paste(paste("Chr", i.chr, sep = ""), i.pos,
                             sep = "@")
    }
    dimnames(result) <- dnames
    class(result) <- "scanqtl"
    attr(result, "method") <- "hk"
    attr(result, "formula") <- qtl:::deparseQTLformula(formula)
    result
}



