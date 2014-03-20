library(qtl)
library(qtl)
library(fda)
source("myfun20130116.R")
library(fields)

# load first replicate of Spalding data
spal <- read.cross("csv", "", "input_rev.csv", genotypes=c("A","B"), na.strings="*")
spal <- convert2riself(spal) # no fcn exist
spal <- jittermap(spal)
#spal <- fill.geno(spal, method=c("imp"))
#spal <- calc.genoprob(spal, step=1)

# drop the two inbred strains
spal <- spal[, c(-(nind(spal)-1), -nind(spal))]

# jitter positions for markers on top of one another
spal <- jittermap(spal)

# grab genetic map
map <- pull.map(spal)

# QTL locations
qtl <- rbind(c(1, 61),
             c(3, 76),
             c(4, 40))

qtlnames <- paste0("QTL", 1:nrow(qtl))

# add QTL positions to the map
mapwqtl <- map
for(i in 1:nrow(qtl)) {
  chr <- qtl[i,1]
  pos <- qtl[i,2]
  mapwqtl[[chr]] <- c(mapwqtl[[chr]], pos)

  # add name
  names(mapwqtl[[chr]])[length(mapwqtl[[chr]])] <- qtlnames[i]

  # sort
  mapwqtl[[chr]] <- sort(mapwqtl[[chr]])
}

# simulate RIL


genD <- function(sampsize) {
    ril <- sim.cross(mapwqtl, n.ind=sampsize, type="riself")

                                        # pull out genotypes for QTL
    qtlgeno <- pull.geno(ril)[, qtlnames]

                                        # drop QTL from cross object
    ril <- drop.markers(ril, qtlnames)

# same pattern of missing data as in spalding data set
    for(i in 1:nchr(ril))
        ril$geno[[i]]$data[is.na(spal$geno[[i]]$data)] <- NA

#dim(ril)
#spal

    m.beta <- c( -0.2677979, -264.5592720,  228.2765707,  -59.4271070)
    covv <- cbind( c(58.99011,  -177.7696,   185.106,   -45.43739),
                  c(-177.76965,  3848.7050, -7274.829,  3595.37470),
                  c(185.10603, -7274.8290, 16897.558, -9702.32353),
                  c(-45.43739,  3595.3747, -9702.324,  6096.71405) )


    Q1 <- c(0.2125015,  8.7758669,  1.2101953, -8.7782531)
    Q2 <- c(-1.909488,  3.395117, -4.683843,  2.680493)
    Q3 <- c(2.481852,   8.988952, -23.307630,  12.737423)


    poly <- function(beta,tt)  # legendre polynomials
        beta[1] + beta[2]*tt + beta[3]*tt^2 + beta[4]*tt^3


    AA <- (qtlgeno - 1.5 )*2


    ttt <- 1:241

    xi <- rnormMulti(sampsize,covv)
################# To control
    rnormMulti <-
        function(n, V)
            matrix(rnorm(ncol(V)*n), ncol=ncol(V)) %*% chol(V)


    Dat <- NULL;
    for( i in 1:sampsize) {
        Dat <- rbind(Dat,
                     m.beta + Q1 * AA[i,1] + Q2 * AA[i,2] + Q3*AA[i,3] + xi[i,] )
    }

    pheno1 <- NULL;
    for(i in 1:sampsize) {
        pheno1 <- rbind(pheno1, poly(Dat[i,],ttt) + rnorm(241, sd = 1) )
    }

    ril$pheno <- data.frame(pheno1)

    ril
}

samples <- genD(162)

simspal <- samples




save(simspal,file = "../data/simspal.RData")


# load first replicate of Spalding data

library(funqtl)
data(exd)
exd <- calc.genoprob(exd, step=1)
out1 <- scanoneF(exd)
perm1 <- scanoneF(exd, n.perm=50)



perm2 <- scanone(exd, pheno.col=1:2, n.perm=100)
summary(perm1)

summary(out1)



class(perm2o)

quantile(perm1[,1], .95)
quantile(perm1[,2], .95)

percentile
str(perm1)

plot(out1)


?scanone



exd
out1 <- scanone(exd, pheno.col=1:10 , method="hk")

out1



data(simspal)
simspal <- calc.genoprob(simspal, step=1)
out1 <- scanone(simspal, pheno.col=1:241 , method="hk")


dim(out1)
z1 <- t(out1[1:661,3:243])
plot2lod(simspal,z1, col=heat.colors(100)[100:1], main="The LOD image of the RIL-1")



data(simspal)
simspal <- calc.genoprob(simspal, step=1)
out1 <- scanone(simspal, pheno.col=1:241 , method="hk")


dim(out1)
z1 <- t(out1[1:234,3:243])
z1 <- t(out1[1:661,3:243])

plot2lod(simspal,out1, col=heat.colors(100)[100:1], main="The LOD image of the RIL-1")


cross <- simspal
z2 <- z1

dim(z2)

out1

table(out2[,2])





  out <- scanone(spal, pheno.col = 1:241, method="hk")


phe <- as.matrix(spal$pheno)
  eff <- NULL
  for(i in 1:nchr(spal)) {
    pr <- spal$geno[[i]]$prob[,,2]
    eff <- rbind(eff, t(apply(pr, 2, function(a,b) lm(b~a)$coef[2], phe)))
  }
  save(out, eff, file=file)



data(simspal)
simspal <- calc.genoprob(simspal, step=1)
out <- scanone(simspal, pheno.col=1:241 , method="hk")
eff <- geteffects(simspal, pheno.cols=1:241)
plotlod(out, eff, gap=15)






dim(eff)
dim(templod)
  templod <- as.matrix(out[,-(1:2)])



dim(eff)

output <- out
effects <- eff


str(as.matrix(effects[[1]]))


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










plot(1:241,simspal$pheno[3,], type = "l", ylim = c(-190,0))
plot(1:241,simspal$pheno[3,], type = "l")


dim(z1)

z1[1:10,1:10]



par(mfrow=c(3,1))
par(las=1)
plot2lod(ril,z1,mval=6.3, col=heat.colors(100)[100:1], main="The LOD image of the RIL-1")
plot2lod(ril2,z2,mval=6.3, col=heat.colors(100)[100:1], main="The LOD image of the RIL-2")
plot2lod(nil3,z3,mval=6.3, col=heat.colors(100)[100:1], main="The LOD image of the NIL data set")


qtl1.c <- makeqtl(ril, chr = c(1, 3, 4),
               pos = c(64, 17, 40.3), what = "prob")

qtl2.c <- makeqtl(ril2, chr = c(1, 3, 4),
               pos = c(64, 17, 40.3), what = "prob")

qtl3.c <- makeqtl(nil3, chr = c(1, 3, 4),
               pos = c(64, 17, 40.3), what = "prob")

thisqtl1.c <- refineqtlF(ril, pheno.cols = pheno.cols,  qtl = qtl1.c,
                      formula = y~Q1 + Q2 + Q3, method = "hk")
thisqtl2.c <- refineqtlF(ril2, pheno.cols = pheno.cols,  qtl = qtl2.c,
                      formula = y~Q1 + Q2 + Q3, method = "hk")
thisqtl3.c <- refineqtlF(nil3, pheno.cols = pheno.cols,  qtl = qtl3.c,
                      formula = y~Q1 + Q2 + Q3, method = "hk")

And ``profileLodMatfn2'' function just make a profile lod lists that
need ``plotlodmatlist2'' function in plotting profile lod scores.
<<>>=
lodmat1.o <- profileLodMatfn2(ril,qtl =  qtl1.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)

lodmat2.o <- profileLodMatfn2(ril2,qtl = qtl2.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)

lodmat3.o <- profileLodMatfn2(nil3,qtl = qtl3.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)

lodmat1.c <- profileLodMatfn2(ril,qtl =  thisqtl1.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)

lodmat2.c <- profileLodMatfn2(ril2,qtl =  thisqtl2.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)

lodmat3.c <- profileLodMatfn2(nil3,qtl =  thisqtl3.c, pheno.cols =1:ncol(ril$pheno),
                             formula = y~Q1 + Q2 + Q3, method = "hk",
                             verbose = F)


par(mfrow=c(3,1))
plotlodmatlist2(lodmat1.o, mval = 8, col=heat.colors(100)[100:1], main="The Profile LOD image of RIL-1")
plotlodmatlist2(lodmat2.o, mval = 8, col=heat.colors(100)[100:1], main="The Profile LOD image of RIL-2")
plotlodmatlist2(lodmat3.o, mval = 5, col=heat.colors(100)[100:1], main="The Profile LOD image of NIL")




scanoneF <-
function(cross, pheno.cols, n.perm, ...) {

    n = nind(cross)

    if (missing(pheno.cols))
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
        for(rep in 1:n.perm)   {
            temp$pheno <- pheno[sample(n),]
            temp <- calc.genoprob(temp, step=0)

            out <- scanone(temp, pheno.col = pheno.cols, ...)
            SLOD <- rowMeans(out[,-(1:2)])
            MLOD <- apply(out[,-(1:2)], 1, max)

            Slods <- c(Slods, max(SLOD) )
            Mlods <- c(Mlods, max(MLOD) )

        }
        return( rbind(Slods,Mlods) )

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
