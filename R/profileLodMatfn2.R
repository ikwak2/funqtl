profileLodMatfn2 <- function (cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL,
    formula, method = c("imp", "hk"), model = c("normal", "binary"),
    verbose = TRUE, tol = 1e-04, maxit.fitqtl = 1000 ) {

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
            else qtl <- makeqtl(cross, chr = chr, pos = pos, what = "prob")
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




    if (!missing(pheno.cols))
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

    formula <- qtl::checkformula(formula, qtl$altname, colnames(covar))


  # identify which QTL are in the model formula
    tovary <- sort(qtl::parseformula(formula, qtl$altname, colnames(covar))$idx.qtl)
    if(length(tovary) != qtl$n.qtl)
        reducedqtl <- qtl::dropfromqtl(qtl, index=(1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl

  # if a QTL is missing from the formula, we need to revise the formula, moving
  # everything over, for use in scanqtl
    if(any(1:length(tovary) != tovary)) {
        tempform <- strsplit(qtl::deparseQTLformula(formula), " *~ *")[[1]][2]
        terms <- strsplit(tempform, " *\\+ *")[[1]]
        for(j in seq(along=terms)) {
            if(length(grep(":", terms[j])) > 0) { # interaction
                temp <- strsplit(terms[j], " *: *")[[1]]

                for(k in seq(along=temp)) {
                    g <- grep("^[Qq][0-9]+$", temp[k])
                    if(length(g) > 0) {
                        num <- as.numeric(substr(temp[k], 2, nchar(temp[k])))
                        temp[k] <- paste("Q", which(tovary == num), sep="")
                    }
                }
                terms[j] <- paste(temp, collapse=":")
            }
            else {
                g <- grep("^[Qq][0-9]+$", terms[j])
                if(length(g) > 0) {
                    num <- as.numeric(substr(terms[j], 2, nchar(terms[j])))
                    terms[j] <- paste("Q", which(tovary == num), sep="")
                }
            }
        }
        formula <- as.formula(paste("y ~", paste(terms, collapse=" + ")))
    }

    curpos <- pos[tovary]
    chrnam <- chr[tovary]

    lc <- length(chrnam)

    lastout <- vector("list", length(curpos))
    names(lastout) <- qtl$name[tovary]

    sexpgm <- getsex(cross)
    cross.attr <- attributes(cross)


    #######
    ######### outline of result form
    ###
    outout <- NULL;
    for(ii in 1:lc) {
        if(!(chrnam[ii] %in% names(outout))) {
            outout[[chrnam[ii]]]$nqtl <- 1
            outout[[chrnam[ii]]]$qtls <- curpos[ii]
            outout[[chrnam[ii]]]$index <- ii
        } else {
            outout[[chrnam[ii]]]$nqtl <- outout[[chrnam[ii]]]$nqtl + 1
            outout[[chrnam[ii]]]$qtls <- c(outout[[chrnam[ii]]]$qtls, curpos[ii])
            outout[[chrnam[ii]]]$index <- c(outout[[chrnam[ii]]]$index, ii)
        }
    }

    for(ii in names(outout) ) {
        outout[[ii]]$index <- outout[[ii]]$index[ order(outout[[ii]]$qtls) ]
    }

#######
#####



    basefit <- NULL
    basefitlod <- NULL

    for(phv in pheno.cols) {
        basefit[[phv]] <- qtl::fitqtlengine(pheno = pheno[,pheno.cols[phv]], qtl = reducedqtl,
                                covar = covar, formula = formula, method = method,
                                model = model, dropone = TRUE, get.ests = FALSE,
                                run.checks = FALSE, cross.attr = cross.attr,
                                sexpgm = sexpgm, tol = tol, maxit = maxit.fitqtl)
        basefitlod <- c( basefitlod, basefit[[phv]]$result.full[1,4] )

    }

    for (j in 1:lc) { # j=1
        otherchr <- chrnam[-j]
        otherpos <- curpos[-j]
        thispos <- as.list(curpos)
        if (any(otherchr == chrnam[j])) {
            linkedpos <- otherpos[otherchr == chr[j]]
            if (any(linkedpos < curpos[j]))
                low <- max(linkedpos[linkedpos < curpos[j]])
            else low <- -Inf
            if (any(linkedpos > curpos[j]))
                high <- min(linkedpos[linkedpos > curpos[j]])
            else high <- Inf
            thispos[[j]] <- c(low, high)
        }  else
        thispos[[j]] <- c(-Inf, Inf)

        for( tt in 1:length(pheno.cols) ) { # tt = 1
            out <- scanqtl(cross = cross, pheno.col = tt,
                           chr = chrnam, pos = thispos, covar = covar, formula = formula,
                           method = method, model = model, incl.markers =TRUE,
                           verbose = scanqtl.verbose, tol = tol, maxit = maxit.fitqtl)

            dropresult <- basefit[[tt]]$result.drop

            if(is.null(dropresult)) {
                if(length(lastout)==1) {
                    dropresult <- rbind(c(NA,NA, basefit$result.full[1,4]))
                    rownames(dropresult) <- names(lastout)
                }
                else
                    stop("There's a problem: need dropresult, but didn't obtain one.")
            }

            rn <- rownames(dropresult)
            qn <- names(lastout)

            if(tt == 1) {
#                names(lastout) <- qtl$name[tovary]
                pos <- as.numeric(matrix(unlist(strsplit(names(out), "@")),
                                         byrow=TRUE,ncol=2)[,2])
                chr <- as.numeric(rep(qtl$chr[tovary][j], length(pos)))

                lastout[[j]] <- cbind(chr, pos,
                        out - (basefit[[tt]][[1]][1,4] - dropresult[rn==qn[j],3])  )
            } else {
                lastout[[j]] <- cbind(lastout[[j]],
                        out - (basefit[[tt]][[1]][1,4] - dropresult[rn==qn[j],3])  )
            }
        }
    }


    #############
    ###########

    for(ii in names(outout)) {
        stp = 1
        for(ij in outout[[ii]]$index) {
            if (stp == 1) {
                outout[[ii]]$out <- lastout[[ij]]
                stp = 2
            } else {
                stposit = 10000000;
                for( iz in nrow(outout[[ii]]$out):1 ) {
                    if ( outout[[ii]]$out[iz,2] == lastout[[ij]][1,2] ) {
                        stposit = iz;
                        break;
                    }
                }
                lst = 1;
                for(iz in stposit:nrow(outout[[ii]]$out))  {
                    outout[[ii]]$out[iz,] <- pmax(outout[[ii]]$out[iz,], lastout[[ij]][lst,])
                    lst = lst + 1;
                }
                if ( lst <= nrow(lastout[[ij]]) ) {
                    outout[[ii]]$out <- rbind(outout[[ii]]$out,
                                              lastout[[ij]][lst:nrow(lastout[[ij]]),] )
                }
            }
        }
    }

###########
    #########

########


########

    for( i in seq(along=lastout)) {
        colnames(lastout[[i]]) <- c("pos", "chr", colnames(cross$pheno))
        lastout[[i]] <- as.data.frame(lastout[[i]])
    }


    # make the profiles scanone objects
    for(i in seq(along=lastout)) { #i =1
        class(lastout[[i]]) <- c("scanmult", "data.frame")
        thechr <- qtl$chr[i]
        if(method=="imp")
            detailedmap <- attr(cross$geno[[thechr]]$draws,"map")
        else
            detailedmap <- attr(cross$geno[[thechr]]$prob,"map")

        if(is.matrix(detailedmap)) detailedmap <- detailedmap[1,]

        r <- range(lastout[[i]][,2])+c(-1e-5, 1e-5)
        rn <- names(detailedmap)[detailedmap>=r[1] & detailedmap<=r[2]]
        o <- grep("^loc-*[0-9]+",rn)
        if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
            rn[o] <- paste("c",thechr,".",rn[o],sep="")
                                        #      if(length(rn) != nrow(lastout[[i]])) return(list(lastout[[i]], rn, detailedmap))
        if(length(rn) == nrow(lastout[[i]])) rownames(lastout[[i]]) <- rn
    }

#    attr(qtl, "lodprofileM") <- lastout
    attr(qtl, "lodprofileM2") <- outout
    class(outout) <- c("lodprofileM2","list")
    outout

}


