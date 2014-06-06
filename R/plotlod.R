#' Plot signed or plain LOD image.
#'
#' Plot signed or plain LOD image. (taking -LOD when the QTL effect
#' is negative), where (i,j) element is for phenotype i at marker j.
#'
#'
#' @param output An object of class \code{"scanone"} as produced by
#' \code{\link[qtl]{scanone}}.
#' @param effects The sign information whether the QTL having AA do negative
#' effect.  Get this by using \code{\link{geteffects}}.
#' @param y Positions of phenotypes in image (e.g., the times).
#' @param ylab y-axis label
#' @param gap The gap between chromosomes.
#' @param ncolors The number of colors between blue and red.
#' @param horizontal If TRUE, transpose the x and y axes
#' @param \dots More graphical components, passed \code{\link[fields]{image.plot}}.
#' @return None.
#' @author Karl W Broman, Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso \code{\link{geteffects}}
#' @keywords hplot
#' @import qtl fields
#' @export
#' @examples
#' data(simspal)
#' simspal <- calc.genoprob(simspal)
#' phe <- 1:nphe(simspal)
#' \dontshow{phe <- seq(1, nphe(simspal), by=60)}
#' out <- scanone(simspal, pheno.col=phe, method="hk")
#' eff <- geteffects(simspal, pheno.cols=phe)
#' nam <- phenames(simspal)
#' y <- as.numeric(substr(nam, 2, nchar(nam)))/60
#' \dontshow{y <- y[phe]}
#' plotlod(out, eff, y,  gap=15)
#' plotlod(out, y=y, gap=15, horizontal = TRUE)
#'
plotlod <- function(output, effects, y, ylab="Time", gap=25,
                    ncolors=251, horizontal=FALSE, ...)
{

    if(missing(y)) {
        y <- 1:(ncol(output)-2)
    }

    templod <- as.matrix(output[,-(1:2)])
    maxlod <- max(templod, na.rm=TRUE)
    zlim <- c(0, maxlod)
    val <- sqrt(seq(0, 1, len=ncolors))
    col <- rgb(1, rev(val), rev(val))
    if(!missing(effects)) {
        if(!all(dim(effects) == dim(templod)))
            stop("dim(effects) doesn't conform to dim(output)")
        templod[effects < 0] <- templod[effects<0] * -1
        zlim <- c(-maxlod, maxlod)
        ncolors=(ncolors-1)/2+1
        val <- sqrt(seq(0, 1, len=ncolors))
        col <- c(rgb(val, val, 1), rgb(1, rev(val), rev(val))[-1])
    }

    uchr <- unique(output[,1])
    pos <- NULL
    lod <- NULL
    chr <- vector("list", length(uchr))
    names(chr) <- uchr
    off.end <- 0.5 # spacing off the ends of the chromosomes
    for(i in seq(along=uchr)) {
        temppos <- output[output[,1]==uchr[i],2]
        temppos <- temppos - min(temppos)
        temppos <- rowMeans(cbind(c(temppos[1]-off.end, temppos),
                                  c(temppos, max(temppos)+off.end)))
        if(is.null(pos)) {
            pos <- temppos
            lod <- templod[output[,1]==uchr[i],]
        }
        else {
            temppos <- max(pos)+gap + temppos
            pos <- c(pos, temppos)
            lod <- rbind(lod, rep(0, ncol(lod)), templod[output[,1]==uchr[i],])
        }
        chr[[i]] <- c(min(temppos), max(temppos))
        chr[[i]] <- c(chr[[i]], mean(chr[[i]]))
    }




    if(horizontal == FALSE) {

        ## add more space between regend and plot
        pos <- c(pos, max(pos) + gap)
        lod <- rbind(lod, rep(0, ncol(lod)))

        ## plot
        image.plot(pos, y, lod, xaxt="n", ylab=ylab, xlab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n")
        title(xlab="Chromosome", mgp=c(2, 0, 0))
        u <- par("usr")
        yd <- 0.04
        for(i in seq(along=chr)) {
            rect(chr[[i]][1]-.5, u[3], chr[[i]][2]+.5, u[3]-diff(u[3:4])*yd, col="gray40", xpd=TRUE)
            text(chr[[i]][3], u[3]-diff(u[3:4])*yd/2, uchr[i], col="white", xpd=TRUE)
            rect(chr[[i]][1]-.5, u[3], chr[[i]][2]+.5, u[4], xpd=TRUE)
        }

    } else {

        image.plot(y, pos, t(lod), yaxt="n", xlab="Time (hours)", ylab="",
              zlim=zlim, col=col, mgp=c(2.6, 1, 0), bty="n")
        title(ylab="Chromosome", mgp=c(2, 0, 0))
        u <- par("usr")
        xd <- 0.04
        for(i in seq(along=chr)) {
            rect(u[1] - diff(u[1:2])*xd, chr[[i]][1]-.5 , u[1], chr[[i]][2]+.5, col="gray40", xpd=TRUE)
            text(u[1]-diff(u[1:2])*xd/2, chr[[i]][3], uchr[i], col="white", xpd=TRUE)
            rect(u[1], chr[[i]][1]-.5 , u[2], chr[[i]][2]+.5, xpd=TRUE)
        }
    }
}

