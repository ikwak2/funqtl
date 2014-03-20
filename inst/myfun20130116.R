library(fields)
library(qtl)




#spal <- calc.genoprob(spal, step=0)
#out.sp <- scanone(spal, pheno.col=1:164, method="hk")
#spal2 <- calc.genoprob(spal2,step=0)
#out.sp2 <- scanone(spal2, pheno.col=1:164, method="hk")

#dim(out.sp)

#out.sp[1:10,1:10]

#z1 <- t(out.sp[,-(1:2)])
#z2 <- t(out.sp2[,-(1:2)])

#plot2lod(spal,z2,mval=6.3, main="The LOD image of the second data set")

#plot2lod(spal,z1,mval=6.3, main="The LOD image of the first data set")

#z1 <- t(LODS1)
#z2 <- t(LODS2)



# plot image of LOD curves
# if effects provided, used signed LOD scores

geteffects <- function(cross,pheno.cols) {
    if(missing(pheno.cols))
        pheno.cols=1:nphe(cross)

    out <- scanone(cross, pheno.col = pheno.cols, method="hk")
    phe <- as.matrix(cross$pheno)
    eff <- NULL
    for(i in 1:nchr(cross)) {
        pr <- cross$geno[[i]]$prob[,,2]
        eff <- rbind(eff, t(apply(pr, 2, function(a,b) lm(b~a)$coef[2,], phe)))
    }
    list(out,eff)
}

plotlod <-
function(output, effects, y, gap=25, off.end=0.5, ncol=251, ...)
{

  mar <- par("mar")
  on.exit(par(mar=mar))
  par(mar=c(3.1, 3.6, 0.6, 0.6), las=1)

  if(missing(y)) { # assume colnames in output are T#
    nam <- names(output)[-(1:2)]
    y <- as.numeric(substr(nam, 2, nchar(nam)))/60
  }

  templod <- as.matrix(output[,-(1:2)])
  maxlod <- max(templod, na.rm=TRUE)
  zlim <- c(0, maxlod)
  val <- sqrt(seq(0, 1, len=ncol))
  col <- rgb(1, rev(val), rev(val))
  if(!missing(effects)) {
    if(!all(dim(effects) == dim(templod)))
      stop("dim(effects) doesn't conform to dim(output)")
    templod[effects < 0] <- templod[effects<0] * -1
    zlim <- c(-maxlod, maxlod)
    ncol=(ncol-1)/2+1
    val <- sqrt(seq(0, 1, len=ncol))
    col <- c(rgb(val, val, 1), rgb(1, rev(val), rev(val))[-1])
  }

  uchr <- unique(output[,1])
  pos <- NULL
  lod <- NULL
  chr <- vector("list", length(uchr))
  names(chr) <- uchr
  for(i in seq(along=uchr)) {
    temppos <- output[output[,1]==i,2]
    temppos <- temppos - min(temppos)
    temppos <- rowMeans(cbind(c(temppos[1]-off.end, temppos),
                              c(temppos, max(temppos)+off.end)))
    if(is.null(pos)) {
      pos <- temppos
      lod <- templod[output[,1]==i,]
    }
    else {
      temppos <- max(pos)+gap + temppos
      pos <- c(pos, temppos)
      lod <- rbind(lod, rep(0, ncol(lod)), templod[output[,1]==i,])
    }
    chr[[i]] <- c(min(temppos), max(temppos))
    chr[[i]] <- c(chr[[i]], mean(chr[[i]]))
  }

  layout(cbind(1, 2), width=c(5, 1))
  image(pos, y, lod, xaxt="n", ylab="Time (hours)", xlab="",
        zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n")
  title(xlab="Chromosome", mgp=c(2, 0, 0))
  u <- par("usr")
  yd <- 0.04
  for(i in seq(along=chr)) {
    rect(chr[[i]][1]-0.5, u[3], chr[[i]][2]+0.5, u[3]-diff(u[3:4])*yd, col="gray40", xpd=TRUE)
    text(chr[[i]][3], u[3]-diff(u[3:4])*yd/2, uchr[i], col="white", xpd=TRUE)
    rect(chr[[i]][1]-0.5, u[3], chr[[i]][2]+0.5, u[4], xpd=TRUE)
  }

  par(mar=c(8.1, 0.6, 4.6, 4.1))
  lodscale <- seq(zlim[1], zlim[2], length=length(col))
  image(0, lodscale, rbind(lodscale), xaxt="n", xlab="", yaxt="n",
        col=col, zlim=zlim)

  axis(side=4, at=pretty(lodscale), lab=abs(pretty(lodscale)))
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], xpd=TRUE)
  text(u[2]+diff(u[1:2])*4, maxlod/2, srt=90, "Cvi > Ler", xpd=TRUE)
  text(u[2]+diff(u[1:2])*4, -maxlod/2, srt=90, "Ler > Cvi", xpd=TRUE)
  text(u[2]+diff(u[1:2])*1.7, u[4]+diff(u[3:4])*0.05, "LOD", xpd=TRUE)

}









plot2lod <- function(cross, z2, ylab="QTL position", xlab="Time", mval=0, col=heat.colors(100)[100:1], cex.n = 1, ...) {
    nch <- nchr(cross)
    nchs <- lapply(cross$geno, function(t) ncol(t$data))
    end <- NULL
    start <- NULL
    num = 0
    z3 <- NULL;

    for(i in 1:nch) { # i=2
        z3 <- cbind(z3, cbind(z2[,((num+1):(num+nchs[[i]]))], rep(mval,nrow(z2)), rep(mval,nrow(z2))))

        num = num + nchs[[i]]
        end <- c(end, num + (i-1)*2 )
        start <- c(start, num - nchs[[i]] + 1 + (i-1)*2 )

    }
    midpt <- (start + end)/2

    x <- 1:(dim(z3)[1])
    y <- 1:(dim(z3)[2])

    par(mar=c(5.1,6.1,2.1,2.1))
    image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab, col=col, ...)

    last <- end[nch]
    start <- c(start, last + 2)
    for(i in 1:nch)
        rect(0,end[i], nrow(z2)+1,start[i+1], col="white")

    rect(0,end[i], nrow(z2)+1,start[i+1]+1, col="white")


    u <- par("usr") # plot ranges [left,right, bottom,top]

    width <- 0.01*diff(u[1:2])
    text(u[1]-5*width, midpt, 1:nch, xpd=TRUE, cex = cex.n)

    for(i in seq(along=start)) {
        if(i %% 2)
            rect(u[1]-width, start[i], u[1]-width*2, end[i], col="gray30", xpd=TRUE)
        else
            rect(u[1]-width*2, start[i], u[1]-width*3, end[i], col="gray30", xpd=TRUE)
    }
}




#cross <- spal2; z1 <- z2 ; mval = 8; lodmatlist <- lodmat2.2; qtl <- thisqtl2.2;

#lodmat2.2.2 <- profileLodMatfn2(spal2,qtl =  thisqtl2.2, pheno.cols =1:ncol(spal$pheno),
#                             formula = y~Q1 + Q2 + Q3 + Q4, method = "hk",
#                             verbose = F)



#lodmatlist <- lodmat2.2.2



#plotlodmatlist2(spal, lodmat2.2.2, mval = 8)


plotlodmatlist2 <- function(lodmatlist, ylab="QTL position", xlab="Time", mval=0, ...) {


    nlst <- length(lodmatlist)
    nlsts <- lapply(lodmatlist, function(t) nrow(t$out))
    chrs <- lapply(lodmatlist, function(t) t$out[1,1])

    end <- NULL
    start <- NULL
    num = 0
    nt <- ncol(lodmatlist[[1]]$out) - 2


    z3 <- NULL;
    qtlpos <- NULL;
    for(i in 1:nlst) { # i=2

        lods <- as.matrix(lodmatlist[[i]]$out[,-c(1,2)])
        z3 <- rbind(z3, rbind(lods, rep(mval,nt), rep(mval,nt)))

        num = num + nlsts[[i]]
        end <- c(end, num + (i-1)*2 )
        start <- c(start, num - nlsts[[i]] + 1 + (i-1)*2 )



        if( lodmatlist[[i]]$nqtl != 0 ) {

            dista = rep(10000, lodmatlist[[i]]$nqtl)
            qtls <- rep(1, lodmatlist[[i]]$nqtl)

            for(j in 1:nlsts[[i]]) {

                for(k in 1:lodmatlist[[i]]$nqtl) {
                    if( abs(lodmatlist[[i]]$out[j,2] - lodmatlist[[i]]$qtls[k]) < dista[k] ) {
                        qtls[k] <- j
                        dista[k] <- abs(lodmatlist[[i]]$out[j,2] - lodmatlist[[i]]$qtls[k])
                    }
                }
            }

            qtlpos <- c(qtlpos, qtls + start[i] - 1)
        }
    }

    midpt <- (start + end)/2

    z3 <- t(z3)
    x <- 1:(dim(z3)[1])
    y <- 1:(dim(z3)[2])

    par(mar=c(5.1,6.1,2.1,2.1))
    image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab, ...)

#        image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab)


    last <- end[nlst]
    start <- c(start, last + 2)
    for(i in 1:nlst)
        rect(0,end[i], nt+1,start[i+1], col="white")

    rect(0,end[i], nt+1,start[i+1]+1, col="white")


    u <- par("usr") # plot ranges [left,right, bottom,top]

    width <- 0.01*diff(u[1:2])
    text(u[1]-5*width, midpt, unlist(chrs)   , xpd=TRUE)

    for(i in seq(along=start)) {
        if(i %% 2)
            rect(u[1]-width, start[i], u[1]-width*2, end[i], col="gray30", xpd=TRUE)
        else
            rect(u[1]-width*2, start[i], u[1]-width*3, end[i], col="gray30", xpd=TRUE)
    }

    for(i in qtlpos) {

        arrows( u[2]+width*2 , i  , u[2] + width/2, i, length = .02, xpd=TRUE )
    }

}


extmat <- function(matlist, schr, pheno.cols, ...) {
    for (ch in schr) {
        matlist[[as.character(ch)]]$nqtl = 0

        ou <- NULL
        k <- 1
        for (i in pheno.cols) {
            outo <- addqtl(chr = ch, pheno.col = i, ...)
            if(k == 1) {
                nn <- nrow(outo)
                ou <- cbind(rep(ch,nn), outo$pos, outo$lod)
                rownames(ou) <- rownames(outo)
                colnames(ou) <- colnames(outo)
                k = 2
            } else {
                ou <- cbind(ou, outo$lod)
            }
        }
        matlist[[as.character(ch)]]$out <- ou
    }
#    matlist
    chrs <- as.numeric(lapply(matlist, function(t) t$out[1,1]))
    oo <- order(chrs)

    matt <- NULL
    for( i in oo) {
        matt[[chrs[i]]] <- matlist[[i]]
    }
    matt
}



plotlodmatlist <- function(lodmatlist, ylab="QTL position", xlab="Time", mval=0, ...) {


    nlst <- length(lodmatlist)
    nlsts <- lapply(lodmatlist, function(t) nrow(t))
    chrs <- lapply(lodmatlist, function(t) t[1,1])
    end <- NULL
    start <- NULL
    num = 0
    nt <- ncol(lodmatlist[[1]]) - 2


    z3 <- NULL;
    for(i in 1:nlst) { # i=1

        lods <- as.matrix(lodmatlist[[i]][,-c(1,2)])
        z3 <- rbind(z3, rbind(lods, rep(mval,nt), rep(mval,nt)))

        num = num + nlsts[[i]]
        end <- c(end, num + (i-1)*2 )
        start <- c(start, num - nlsts[[i]] + 1 + (i-1)*2 )

    }
    midpt <- (start + end)/2

    z3 <- t(z3)
    x <- 1:(dim(z3)[1])
    y <- 1:(dim(z3)[2])

    par(mar=c(5.1,6.1,2.1,2.1))
    image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab, ...)

#        image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab)


    last <- end[nlst]
    start <- c(start, last + 2)
    for(i in 1:nlst)
        rect(0,end[i], nt+1,start[i+1], col="white")

    rect(0,end[i], nt+1,start[i+1]+1, col="white")


    u <- par("usr") # plot ranges [left,right, bottom,top]

    width <- 0.01*diff(u[1:2])
    text(u[1]-5*width, midpt, unlist(chrs)   , xpd=TRUE)

    for(i in seq(along=start)) {
        if(i %% 2)
            rect(u[1]-width, start[i], u[1]-width*2, end[i], col="gray30", xpd=TRUE)
        else
            rect(u[1]-width*2, start[i], u[1]-width*3, end[i], col="gray30", xpd=TRUE)
    }

}









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

    formula <- qtl:::checkformula(formula, qtl$altname, colnames(covar))


  # identify which QTL are in the model formula
    tovary <- sort(qtl:::parseformula(formula, qtl$altname, colnames(covar))$idx.qtl)
    if(length(tovary) != qtl$n.qtl)
        reducedqtl <- qtl:::dropfromqtl(qtl, index=(1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl

  # if a QTL is missing from the formula, we need to revise the formula, moving
  # everything over, for use in scanqtl
    if(any(1:length(tovary) != tovary)) {
        tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
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
        basefit[[phv]] <- qtl:::fitqtlengine(pheno = pheno[,pheno.cols[phv]], qtl = reducedqtl,
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

}










profileLodMatfn <- function (cross, pheno.cols, qtl, chr, pos, qtl.name, covar = NULL,
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

    formula <- qtl:::checkformula(formula, qtl$altname, colnames(covar))


  # identify which QTL are in the model formula
    tovary <- sort(qtl:::parseformula(formula, qtl$altname, colnames(covar))$idx.qtl)
    if(length(tovary) != qtl$n.qtl)
        reducedqtl <- qtl:::dropfromqtl(qtl, index=(1:qtl$n.qtl)[-tovary])
    else reducedqtl <- qtl

  # if a QTL is missing from the formula, we need to revise the formula, moving
  # everything over, for use in scanqtl
    if(any(1:length(tovary) != tovary)) {
        tempform <- strsplit(qtl:::deparseQTLformula(formula), " *~ *")[[1]][2]
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

    attr(qtl, "lodprofileM") <- lastout

}







addqtlF <- function(cross, pheno.cols, ...) {

    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    LODS <- NULL;
    for(i in pheno.cols ) {
        out <- addqtl(cross, pheno.col = i, ...)
        LODS <- cbind(LODS, out$lod)
    }

    MXy <- max(LODS)
    Slods <- apply(LODS, 1, mean)
    Mlods <- apply(LODS, 1, max)
    out[,3] <- Slods
    out[,4] <- Mlods
    names(out)[3:4] <- c("slod","mlod")

    out[,1:4]
}





scanoneF <- function(cross, pheno.cols=NULL, n.perm, ...) {

    n = nind(cross)

    if (missing(pheno.cols)) {
        pheno.cols = 1:nphe(cross)
    }
    if (!all(pheno.cols %in% 1:nphe(cross))) {
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))
    }
    if (missing(n.perm))
        n.perm <- 0

    if (n.perm > 0 ) {

        temp <- cross
        pheno <- cross$pheno[,pheno.cols]

        Slods <- NULL;
        Mlods <- NULL;
        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]
            temp <- calc.genoprob(temp, step=0)

            out <- scanone(temp, pheno.col = pheno.cols, ...)
            SLOD <- rowMeans(out[,-(1:2)])
            MLOD <- apply(out[,-(1:2)], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )

        }
        return( cbind(Slods,Mlods) )

    } else {

        out <- scanone(cross, pheno.col = pheno.cols, ...)
        SLOD <- rowMeans(out[,-(1:2)])
        MLOD <- apply(out[,-(1:2)], 1, max)

        out[,3] <- SLOD
        out[,4] <- MLOD
        names(out)[3:4] <- c("slod","mlod")

        out[,1:4]
    }
}



scantwoF <- function(cross, pheno.cols, usec=c("slod","mlod"), n.perm, ...) {

    n = nind(cross)
    usec <- match.arg(usec)

    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    if (missing(n.perm))
        n.perm <- 0

    if (n.perm > 0 ) {
        temp <- cross
        pheno <- cross$pheno[,pheno.cols]

        Slods <- NULL;
        Mlods <- NULL;
        Slod <- NULL;
        Mlod <- NULL;
        SlodsH <- NULL;
        SlodsL <- NULL;
        MlodsH <- NULL;
        MlodsL <- NULL;

        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]
            temp <- calc.genoprob(temp, step=0)

            out2 <- scantwo(temp, pheno.col = pheno.cols,  ...)
            out1 <- scanone(temp, pheno.col = pheno.cols,  ...)

            # out3 for slod
            out3 <- out2
            out3$lod <- apply(out2$lod, 1:2, mean)

            # out4 for mlod
            out4 <- out2
            out4$lod <- apply(out2$lod, 1:2, max)

            SLOD <- rowMeans(out1[,-(1:2)])
            MLOD <- apply(out1[,-(1:2)], 1, max)

            Slod <- c(Slod, max(summary(out3)$one))
            Mlod <- c(Mlod, max(summary(out4)$one))
            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )
            SlodsH <- c(SlodsH, max(summary(out3)$lod.int) )
            SlodsL <- c(SlodsL, max(summary(out3)$lod.fv1) )
            MlodsH <- c(MlodsH, max(summary(out4)$lod.int) )
            MlodsL <- c(MlodsL, max(summary(out4)$lod.fv1) )
        }

        return( cbind(Slod,Mlod,Slods,Mlods,SlodsH, SlodsL, MlodsH, MlodsL) )

    } else {
        out <- scantwo(cross, pheno.col = pheno.cols, ...)

        if(usec=="slod") {
            out$lod <- apply(out$lod, 1:2, mean)
        }
        if(usec=="mlod") {
            out$lod <- apply(out$lod, 1:2, max)
        }
    out
    }
}








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







addintF <- function (cross, pheno.cols = pheno.cols, ...) {
    if (!("cross" %in% class(cross)))
        stop("The cross argument must be an object of class \"cross\".")
    if (!("qtl" %in% class(qtl)))
        stop("The qtl argument must be an object of class \"qtl\".")
    if (!is.null(covar) && !is.data.frame(covar)) {
        if (is.matrix(covar) && is.numeric(covar))
            covar <- as.data.frame(covar, stringsAsFactors = TRUE)
        else stop("covar should be a data.frame")
    }
    if (LikePheVector(pheno.col, nind(cross), nphe(cross))) {
        cross$pheno <- cbind(pheno.col, cross$pheno)
        pheno.col <- 1
    }
    if (length(pheno.col) > 1) {
        pheno.col <- pheno.col[1]
        warning("addint can take just one phenotype; only the first will be used")
    }
    if (is.character(pheno.col)) {
        num <- find.pheno(cross, pheno.col)
        if (is.na(num))
            stop("Couldn't identify phenotype \"", pheno.col,
                "\"")
        pheno.col <- num
    }
    if (pheno.col < 1 | pheno.col > nphe(cross))
        stop("pheno.col values should be between 1 and the no. phenotypes")
    pheno <- cross$pheno[, pheno.col]
    if (!is.null(covar) && nrow(covar) != length(pheno))
        stop("nrow(covar) != no. individuals in cross.")
    method <- match.arg(method)
    model <- match.arg(model)
    if (!missing(formula) && is.character(formula))
        formula <- as.formula(formula)
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
    if (!is.null(covar))
        phcovar <- cbind(pheno, covar)
    else phcovar <- as.data.frame(pheno, stringsAsFactors = TRUE)
    if (any(is.na(phcovar))) {
        if (ncol(phcovar) == 1)
            hasmissing <- is.na(phcovar)
        else hasmissing <- apply(phcovar, 1, function(a) any(is.na(a)))
        if (all(hasmissing))
            stop("All individuals are missing phenotypes or covariates.")
        if (any(hasmissing)) {
            warning("Dropping ", sum(hasmissing), " individuals with missing phenotypes.\n")
            pheno <- pheno[!hasmissing]
            qtl$n.ind <- sum(!hasmissing)
            if (method == "imp")
                qtl$geno <- qtl$geno[!hasmissing, , , drop = FALSE]
            else qtl$prob <- lapply(qtl$prob, function(a) a[!hasmissing,
                , drop = FALSE])
            if (!is.null(covar))
                covar <- covar[!hasmissing, , drop = FALSE]
        }
    }
    if (is.null(covar))
        n.covar <- 0
    else n.covar <- ncol(covar)
    n.qtl <- qtl$n.qtl
    if (missing(formula)) {
        tmp.Q <- paste("Q", 1:n.qtl, sep = "")
        formula <- "y~Q1"
        if (n.qtl > 1)
            for (i in 2:n.qtl) formula <- paste(formula, tmp.Q[i],
                sep = "+")
        if (n.covar) {
            tmp.C <- colnames(covar)
            for (i in 1:n.covar) formula <- paste(formula, tmp.C[i],
                sep = "+")
        }
        formula <- as.formula(formula)
    }
    formula <- checkformula(formula, qtl$altname, colnames(covar))
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
    thefit0 <- fitqtlengine(pheno = pheno, qtl = qtl, covar = covar,
        formula = formula, method = method, model = model, dropone = FALSE,
        get.ests = FALSE, run.checks = FALSE, cross.attr = cross.attr,
        sexpgm = sexpgm, tol = tol, maxit = maxit)
    results <- matrix(ncol = 7, nrow = n2test)
    dimnames(results) <- list(int2test.alt, c("df", "Type III SS",
        "LOD", "%var", "F value", "Pvalue(Chi2)", "Pvalue(F)"))
    for (k in seq(along = int2test)) {
        thefit1 <- fitqtlengine(pheno = pheno, qtl = qtl, covar = covar,
            formula = as.formula(paste(deparseQTLformula(formula),
                int2test[k], sep = "+")), method = method, model = model,
            dropone = FALSE, get.ests = FALSE, run.checks = FALSE,
            cross.attr = cross.attr, sexpgm = sexpgm, tol = tol,
            maxit = maxit)
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
    attr(results, "method") <- method
    attr(results, "model") <- model
    attr(results, "formula") <- deparseQTLformula(formula)
    if (simple)
        pvalues <- FALSE
    attr(results, "pvalues") <- pvalues
    attr(results, "simple") <- simple
    results
}










stepwiseqtlF <- function (cross, chr, pheno.cols, qtl, usec=c("slod","mlod"), formula, max.qtl = 10,
    covar = NULL, method = c("imp", "hk"), model = c("normal",
        "binary"), incl.markers = TRUE, refine.locations = TRUE,
    additive.only = FALSE, penalties, keeplodprofile = FALSE,
    keeptrace = FALSE, verbose = TRUE, tol = 1e-04, maxit = 1000)
{

    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    #
    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))


    pheno <- cross$pheno[,pheno.cols]

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
        if (keeplodprofile) {
            if (verbose)
                cat(" ---One last pass through refineqtl\n")
            qtl <- refineqtl(cross, pheno.col = pheno.col, qtl = qtl,
                covar = covar, formula = formula, method = method,
                verbose = verbose.scan, incl.markers = incl.markers,
                             keeplodprofile = TRUE)
        }
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














#scanqtlfn <- function (cross, pheno.cols, chr, pos, covar = NULL, formula,
#    method = c("imp", "hk"), model = c("normal", "binary"), incl.markers = FALSE,
#    verbose = TRUE, tol = 1e-04, maxit = 1000)

scanqtlfn <- function (cross, pheno.cols, chr, pos, covar = NULL, formula,
    method = c("imp", "hk"), model = c("normal", "binary"), incl.markers = FALSE,
    verbose = TRUE, tol = 1e-04, maxit = 1000)
{
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if (!is.null(covar) && !is.data.frame(covar)) {
        if (is.matrix(covar) && is.numeric(covar))
            covar <- as.data.frame(covar, stringsAsFactors = TRUE)
        else stop("covar should be a data.frame")
    }
#    if (qtl:::LikePheVector(pheno.col, nind(cross), nphe(cross))) {
#        cross$pheno <- cbind(pheno.col, cross$pheno)
#        pheno.col <- 1
#    }
#    if (length(pheno.col) > 1) {
#        pheno.col <- pheno.col[1]
#        warning("scanqtl can take just one phenotype; only the first will be used")
#    }
#    if (is.character(pheno.col)) {
#        num <- find.pheno(cross, pheno.col)
#        if (is.na(num))
#            stop("Couldn't identify phenotype \"", pheno.col,
#                "\"")
#        pheno.col <- num
#    }

    if (!missing(pheno.cols))
        pheno.cols = 1:nphe(cross)

    #
    if (!all(pheno.cols %in% 1:nphe(cross)))
        stop("pheno.cols should be in a range of 1 to ", nphe(cross))

    #
    pheno <- as.data.frame(cross$pheno[, pheno.cols], stringsAsFactors = TRUE)

    if (!is.null(covar) && nrow(covar) != length(pheno))
        stop("nrow(covar) != no. individuals in cross.")

##########
###########
#    method <- match.arg(method)
#    model <- match.arg(model)


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
#    method <- match.arg(method)
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
        formula.str <- qtl:::deparseQTLformula(formula)
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
    formula <- qtl:::checkformula(formula, paste("Q", 1:length(chr),
        sep = ""), colnames(covar))
    if (!is.null(covar)) {
        theterms <- rownames(attr(terms(formula), "factors"))
        m <- match(colnames(covar), theterms)
        if (all(is.na(m)))
            covar <- NULL
        else covar <- covar[, !is.na(m), drop = FALSE]
    }
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
#
            pheno <- pheno[!hasmissing,]
            if (!is.null(covar))
                covar <- covar[!hasmissing, , drop = FALSE]
        }
    }
    sexpgm <- getsex(cross)
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
            } #i=4
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


        results <- NULL
        results2 <- NULL
        for (ii in 1:length(pheno.cols)) {
            results[[ii]] <- qtl:::fitqtlengine(pheno = pheno[,pheno.cols[ii]], qtl = qtl, covar = covar,
            formula = formula, method = method, model = model,
            dropone = FALSE, get.ests = FALSE, run.checks = FALSE,
            cross.attr = cross.attr, sexpgm = sexpgm, tol = tol,
            maxit = maxit)

            results2 <- c(results2, results[[ii]][[1]][1,4])
        }
        result <- mean(results2)  #result[[1]][1, 4]



        names(result) <- "MLOD"
        class(result) <- "scanqtlfn"
        attr(result, "method") <- method
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

##

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
#
                        qtl.obj$geno[, kk, ] <- qtl:::reviseXdata(type,
                        "full", sexpgm, draws = qtl.obj$geno[,
                          kk, , drop = FALSE], cross.attr = attributes(cross))
                    else {
                      temp <- qtl.obj$prob[[kk]]
                      temp <- array(temp, dim = c(nrow(temp),
                        1, ncol(temp)))
#
                      dimnames(temp) <- list(NULL, "loc", 1:ncol(qtl.obj$prob[[kk]]))
                      qtl.obj$prob[[kk]] <- qtl:::reviseXdata(type,
                        "full", sexpgm, prob = temp, cross.attr = attributes(cross))[,
                        1, ]
                    }
                  }
                  current.pos[kk] <- pos.tmp[kk]
                }
            }
        }


#
        fit <- NULL;
        fitresults <- NULL;
        for(ii in 1:length(pheno.cols)) {
            fit[[ii]] <- qtl:::fitqtlengine(pheno = pheno[,pheno.cols[ii]], qtl = qtl.obj, covar = covar,
                 formula = formula, method = method, model = model,
                 dropone = FALSE, get.ests = FALSE, run.checks = FALSE,
                 cross.attr = cross.attr, sexpgm = sexpgm, tol = tol,
                 maxit = maxit)
            fitresults <- c(fitresults, fit[[ii]][[1]][1,4])
        }

        if (verbose && ((i - 1)%%n.prnt) == 0)
            cat("    ", i, "/", n.loop, "\n")
#
        result[i] <- mean(fitresults)
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
    attr(result, "formula") <- qtl:::deparseQTLformula(formula)
    result
}

gen.data2 <- function(sample.size, cov.fcn, beta.coef,er){
  ##popu.size <- 10000  # population size

  ## simulate genotypes
  mp <- sim.map(rep(100,4), n.mar=6, include.x=F, eq.spacing=T) # simulate map
  md <- c(1,32,0,0)                       # one QTL at 32cM on chrom. 1
  samples <- sim.cross(map=mp, model=md, type='f2', n.ind=sample.size, keep.qtlgeno=T)

  ## retrieve sample genotypes
  ## ind <- sample(popu.size, sample.size)
  ## samples <- subset.cross(cross, ind=ind)
  ## ## subsetting qtlgeno by hand since subset.cross does not do it.
  ## qtlgeno <- samples$qtlgeno[ind]
  ## samples$qtlgeno <- qtlgeno

  ## simulate phenotypes
  ## 1. get the means
  mean.vals <- matrix(0., nrow=3, ncol=len.tt) # only 3 genotypes means only 3 mean curves
  for (i in 1:3){
    mean.vals[i,] <- logisticFun(beta.coef[i,],tt)
  }
  sample.means <- mean.vals[samples$qtlgeno,]
  ## 2. get the noises
  if (cov.fcn == 'autocorr'){           # case (1)
    ee <- rnormAutocor(sample.size, tt, 0.6, 3.0*er) # rho = 0.6, sigma^2 = 3
  }
  else if (cov.fcn == 'equicorr'){      # case (2)
    cov <- matrix(0.5, nrow=len.tt, ncol=len.tt)  # rho = 0.5
    diag(cov) <- 1.0
    cov <- 3.0 * cov                    # sigma^2 = 3.0
    ee <- rnormMulti(sample.size, cov*er)
  }
  else if (cov.fcn == 'structured'){    # case(3)
    ee <- rnormMulti(sample.size, structured.cov*er)
  }
  else{
    stop('Unknown covariance function.')
  }
  Y <- sample.means + ee

  samples$pheno <- Y
  return(samples)
}

compsim2 <- function(iter, N, method, beta.coef,er ) {
    results = NULL;
    resultm = NULL;
    resultss = NULL;
    resultqf = NULL;
    for(rep in 1:iter) {
        D1 <- gen.data2(N, method, beta.coef,er)
        D1 <- calc.genoprob(D1, step=4)
        D1.out <- scanoneF(D1, pheno.cols=1:10, method="hk")

        pheno <- D1$pheno;
        tt <- 1:(ncol(pheno))
        phi5 <-  bs(tt, df=5, intercept = FALSE)
        D1ss.out <- funcScanone(pheno, D1, phi5, crit = "ss")
        D1qf.out <- funcScanone(pheno, D1, phi5, crit = "qf")
        locass <- D1ss.out[which(D1ss.out$lod == max(D1ss.out$lod)), 1:3]
        locaqf <- D1qf.out[which(D1qf.out$lod == max(D1qf.out$lod)), 1:3]

        locas <- D1.out[which(D1.out$slod == max(D1.out$slod)),1:3]
        locam <- D1.out[which(D1.out$mlod == max(D1.out$mlod)),c(1,2,4)]
                                        #        if(loca[1] == 1)
        results <- rbind(results, as.vector(locas))
        resultm <- rbind(resultm, as.vector(locam))
        resultss <- rbind(resultss, as.vector(locass))
        resultqf <- rbind(resultqf, as.vector(locaqf))
    }

    out <- list( slod = results, mlod = resultm, ss = resultss, qf = resultqf)
}


permsim <- function(iter, N, method, beta.coef,er ) {
    results = NULL;
    resultm = NULL;
    resultss = NULL;
    resultqf = NULL;
    for(rep in 1:iter) {
        D1 <- gen.data2(N, method, beta.coef,er)
        D1 <- calc.genoprob(D1, step=4)
        D1.out <- scanoneF(D1, pheno.cols=1:10, method="hk")

        o <- sample(N)

        pheno <- D1$pheno[o,];
        D1$pheno <- pheno
        tt <- 1:(ncol(pheno))
        phi5 <-  bs(tt, df=5, intercept = FALSE)
        D1ss.out <- funcScanone(pheno, D1, phi5, crit = "ss")
        D1qf.out <- funcScanone(pheno, D1, phi5, crit = "qf")
        locass <- D1ss.out[which(D1ss.out$lod == max(D1ss.out$lod)), 1:3]
        locaqf <- D1qf.out[which(D1qf.out$lod == max(D1qf.out$lod)), 1:3]

        locas <- D1.out[which(D1.out$slod == max(D1.out$slod)),1:3]
        locam <- D1.out[which(D1.out$mlod == max(D1.out$mlod)),c(1,2,4)]
                                        #        if(loca[1] == 1)
        results <- rbind(results, as.vector(locas))
        resultm <- rbind(resultm, as.vector(locam))
        resultss <- rbind(resultss, as.vector(locass))
        resultqf <- rbind(resultqf, as.vector(locaqf))
    }

    out <- list( slod = results, mlod = resultm, ss = resultss, qf = resultqf)
}

prtout <- function( res, lodcrit, window ) {
    respos <- res[res[,1] == "1", 2]
    reslod <- res[res[,1] == "1", 3]
    TP = ( abs(respos - 32) < window & reslod > lodcrit )
    FP = ( abs(respos - 32) >= window & reslod > lodcrit )
    FN = ( abs(respos - 32) < window & reslod <= lodcrit )
    TN = ( abs(respos - 32) >= window & reslod <= lodcrit )
    prtout <- list(length =  length(respos),
                   mean = mean(respos),
                   Pmean = mean(respos[reslod > lodcrit]),
                   sd = sd(respos),
                   rmse = sqrt(mean((respos-32)^2)),
                   TP = sum(TP), FP = sum(FP), FN = sum(FN), TN = sum(TN)
                   )
}


