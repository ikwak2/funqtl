scanqtlfn <-
function (cross, pheno.cols, chr, pos, covar = NULL, formula,
          method = c("hk", "imp"), incl.markers = FALSE,
          verbose = TRUE, usec = c("slod", "mlod") )
{
  usec <- match.arg(usec)

  # check inputs
  if (!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")
  if (!is.null(covar) && !is.data.frame(covar)) {
    if (is.matrix(covar) && is.numeric(covar))
      covar <- as.data.frame(covar, stringsAsFactors = TRUE)
    else stop("covar should be a data.frame")
  }

  if (missing(pheno.cols))
    pheno.cols = 1:nphe(cross)

  if (!all(pheno.cols %in% 1:nphe(cross)))
    stop("pheno.cols should be in a range of 1 to ", nphe(cross))

  pheno <- as.data.frame(cross$pheno[, pheno.cols,drop=FALSE], stringsAsFactors = TRUE)

  if (!is.null(covar) && nrow(covar) != nrow(pheno))
    stop("nrow(covar) != no. individuals in cross.")

  # check formula and qtl object
  if (!missing(formula) && is.character(formula))
    formula <- as.formula(formula)
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
  if (method == "imp") {
    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$draws)) &&
        attr(cross$geno[[1]]$draws, "stepwidth") != "fixed") {
      stepwidth.var <- TRUE
      incl.markers <- TRUE
    }
    else stepwidth.var <- FALSE
  }
  else {
    if ("stepwidth" %in% names(attributes(cross$geno[[1]]$prob)) &&
        attr(cross$geno[[1]]$prob, "stepwidth") != "fixed") {
      stepwidth.var <- TRUE
      incl.markers <- TRUE
    }
    else stepwidth.var <- FALSE
  }
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno, class)
  if (length(chr) != length(pos))
    stop("Input chr and pos must have the same length")

  ichr <- match(chr, names(cross$geno))
  if (any(is.na(ichr)))
    stop("There's no chromosome number ", chr[is.na(ichr)],
         " in input cross object")
  n.qtl <- length(chr)
  n.covar <- length(covar)
  if (missing(formula)) {
    tmp.Q <- paste("Q", 1:n.qtl, sep = "")
    formula <- "y~Q1"
    if (n.qtl > 1)
      for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                                          sep = "+")
    if (n.covar) {
      tmp.C <- names(covar)
      for (i in 1:n.covar) formula <- paste(formula, tmp.C[i],
                                            sep = "+")
    }
    formula <- as.formula(formula)
  }
  else {
    formula.str <- qtl::deparseQTLformula(formula)
    for (i in 1:n.qtl) {
      qtl.term <- paste("Q", i, sep = "")
      if (length(grep(qtl.term, formula.str, ignore.case = TRUE)) ==
          0)
        formula.str <- paste(formula.str, qtl.term, sep = "+")
    }
    if (n.covar) {
      for (i in 1:n.covar) {
        covar.term <- names(covar)[i]
        if (length(grep(covar.term, formula.str, ignore.case = TRUE)) ==
            0)
          formula.str <- paste(formula.str, covar.term,
                               sep = "+")
      }
    }
    formula <- as.formula(formula.str)
  }
  formula <- qtl::checkformula(formula, paste("Q", 1:length(chr),
                                              sep = ""), colnames(covar))
  if (!is.null(covar)) {
    theterms <- rownames(attr(terms(formula), "factors"))
    m <- match(colnames(covar), theterms)
    if (all(is.na(m)))
      covar <- NULL
    else covar <- covar[, !is.na(m), drop = FALSE]
  }

  # deal with missing data
  if (!is.null(covar))
    phcovar <- cbind(pheno, covar)
  else phcovar <- pheno
  if (any(is.na(phcovar))) {
    if (ncol(phcovar) == 1)
      hasmissing <- is.na(phcovar)
    else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
    if (all(hasmissing))
      stop("All individuals are missing phenotypes or covariates.")
    if (any(hasmissing)) {
      warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
      cross <- subset(cross, ind = !hasmissing)
      pheno <- pheno[!hasmissing,]
      if (!is.null(covar))
        covar <- covar[!hasmissing, , drop = FALSE]
    }
  }
  sexpgm <- getsex(cross)

  # null LOD
  lod0 <- rep(0, length(pheno.cols))
  if(!is.null(covar)) {
    pheno <- cross$pheno[,pheno.cols,drop=FALSE]
    rss0 <- colSums(lm(as.matrix(pheno) ~ as.matrix(covar))$resid^2, na.rm=TRUE)
    rss00 <- colSums(lm(as.matrix(pheno) ~ 1)$resid^2, na.rm=TRUE)
    lod0 <- nrow(pheno)/2 * log10(rss00/rss0)
  }

  idx.varied <- NULL
  indices <- pos
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
      if (method == "imp") {
        if ("map" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
          map <- attr(cross$geno[[ichr[i]]]$draws, "map")
        else {
          stp <- attr(cross$geno[[ichr[i]]]$draws, "step")
          oe <- attr(cross$geno[[ichr[i]]]$draws, "off.end")
          if ("stepwidth" %in% names(attributes(cross$geno[[ichr[i]]]$draws)))
            stpw <- attr(cross$geno[[ichr[i]]]$draws,
                         "stepwidth")
          else stpw <- "fixed"
          map <- create.map(cross$geno[[ichr[i]]]$map,
                            stp, oe, stpw)
        }
      }
      else {
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
      }
      if (is.matrix(map))
        map <- map[1, ]
      indices[[i]] <- seq(along = map)
      if (method == "imp")
        step <- attr(cross$geno[[ichr[i]]]$draws, "step")
      else step <- attr(cross$geno[[ichr[i]]]$prob, "step")
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
  }
  else {
    if (method == "imp")
      qtl <- makeqtl(cross, chr = chr, pos = unlist(pos),
                     what = "draws")
    else qtl <- makeqtl(cross, chr = chr, pos = unlist(pos),
                        what = "prob")


    fitresults <- rep(NA, length(pheno.cols))
    for (ii in 1:length(pheno.cols)) {
      fit <- qtl::fitqtlengine(pheno = pheno[,ii], qtl = qtl, covar = covar,
                               formula = formula, method = method, model = "normal",
                               dropone = FALSE, get.ests = FALSE, run.checks = FALSE,
                               cross.attr = cross.attr, sexpgm = sexpgm)

      fitresults[ii] <- fit[[1]][1,4]
    }
    result <- ifelse(usec=="slod", mean(fitresults-lod0), max(fitresults-lod0))

    names(result) <- toupper(usec)
    class(result) <- "scanqtlfn"
    attr(result, "method") <- method
    attr(result, "formula") <- qtl::deparseQTLformula(formula)
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
      if (method == "imp")
        qtl.obj <- makeqtl(cross, chr, pos.tmp, what = "draws")
      else qtl.obj <- makeqtl(cross, chr, pos.tmp, what = "prob")
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
          if (method == "imp")
            qtl.obj$geno[, kk, ] <- cross$geno[[ichr[kk]]]$draws[,
                                                                 w, ]
          else qtl.obj$prob[[kk]] <- cross$geno[[ichr[kk]]]$prob[,
                                                                 w, ]
          thew[kk] <- w
          if (chrtype[ichr[kk]] == "X" && (type == "bc" ||
                       type == "f2")) {
            if (method == "imp")

              qtl.obj$geno[, kk, ] <- qtl::reviseXdata(type,
                                                       "full", sexpgm, draws = qtl.obj$geno[,
                                                                         kk, , drop = FALSE], cross.attr = attributes(cross))
            else {
              temp <- qtl.obj$prob[[kk]]
              temp <- array(temp, dim = c(nrow(temp),
                                    1, ncol(temp)))

              dimnames(temp) <- list(NULL, "loc", 1:ncol(qtl.obj$prob[[kk]]))
              qtl.obj$prob[[kk]] <- qtl::reviseXdata(type,
                                                     "full", sexpgm, prob = temp, cross.attr = attributes(cross))[,
                                                                                    1, ]
            }
          }
          current.pos[kk] <- pos.tmp[kk]
        }
      }
    }

    fitresults <- rep(NA, length(pheno.cols))
    for(ii in 1:length(pheno.cols)) {
      fit <- qtl::fitqtlengine(pheno = pheno[,ii], qtl = qtl.obj, covar = covar,
                               formula = formula, method = method, model = "normal",
                               dropone = FALSE, get.ests = FALSE, run.checks = FALSE,
                               cross.attr = cross.attr, sexpgm = sexpgm)
      fitresults[ii] <- fit[[1]][1,4]
    }

    if (verbose && ((i - 1)%%n.prnt) == 0)
      cat("    ", i, "/", n.loop, "\n")
    result[i] <- ifelse(usec=="slod", mean(fitresults-lod0), max(fitresults-lod0))

  }
  dnames <- list(NULL)
  for (i in 1:n.idx.varied) {
    i.chr <- chr[idx.varied[n.idx.varied - i + 1]]
    i.pos <- pos[[idx.varied[n.idx.varied - i + 1]]]
    dnames[[i]] <- paste(paste("Chr", i.chr, sep = ""), i.pos,
                         sep = "@")
  }
  dimnames(result) <- dnames
  class(result) <- "scanqtlfn"
  attr(result, "method") <- method
  attr(result, "formula") <- qtl::deparseQTLformula(formula)
  result
}
