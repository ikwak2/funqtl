plotlod <- function(output, effects, y, gap=25, off.end=0.5, ncol=251, ...)
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
#  u <- par("usr")
#  rect(u[1], u[3], u[2], u[4], xpd=TRUE)
#  text(u[2]+diff(u[1:2])*4, maxlod/2, srt=90, "Cvi > Ler", xpd=TRUE)
#  text(u[2]+diff(u[1:2])*4, -maxlod/2, srt=90, "Ler > Cvi", xpd=TRUE)
#  text(u[2]+diff(u[1:2])*1.7, u[4]+diff(u[3:4])*0.05, "LOD", xpd=TRUE)

}
