
## ----knitr_options, echo=FALSE, results=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.width = 12)


## ----loading_run, echo=FALSE, message=FALSE------------------------------
library(funqtl)
data(simspal)


## ----loading, eval=FALSE-------------------------------------------------
## library(funqtl)
## data(simspal)


## ----summary_simspal-----------------------------------------------------
summary(simspal)


## ----plot_map_and_pheno--------------------------------------------------
par(mfrow=c(1,2))
plotMap(simspal, main="")
plot(1:241, simspal$pheno[160,], type="l", xlab="Time", ylim=c(-120,0),
     ylab="Root Tip Angle (degrees)")
ind <- c(19, 20, 132, 72)
color <- c("blue", "red", "green", "orange")
for(i in seq(along=ind))
  lines(1:241, simspal$pheno[ind[i],], col=color[i])


## ----scanone-------------------------------------------------------------
n.phe <- nphe(simspal)
simspal <- calc.genoprob(simspal, step=0)
out <- scanone(simspal, pheno.col = 1:n.phe, method="hk")


## ----geteffects----------------------------------------------------------
eff <- geteffects(simspal, pheno.cols=1:n.phe)


## ----plotload, fig.height=10, fig.width = 12-----------------------------
plotlod(out, eff, gap=15, main="The LOD image of the simspal data set",
        ylab="Time")


## ----slodmlod------------------------------------------------------------
out1 <- scanoneF(simspal, pheno.cols = 1:241, method="hk")



## ----slodmlodcurve, fig.height=10, fig.width = 12------------------------
par(mfrow=c(2,1))
plot(out1, ylim=c(0,3.5), main="The SLOD curve for IBMset1", bandcol="gray90")
abline(h=2.02, col="red", lty=3)
plot(out1[,c(1,2,4)], ylim=c(0,7), main="The MLOD curve for IBMset1", bandcol="gray90")
abline(h=3.46, col="red", lty=3)

# permutation threshold
# permout <- scanoneF(simspal, pheno.cols=1:241, method = "hk", n.perm=1000)
# summary(permout) # display 5, 10 % threshold of permutation result


## ----stepwiseqtlscan-----------------------------------------------------

#qtlslod <- stepwiseqtlF(simspal, pheno.cols = 1:241, max.qtl = 6,
#                       usec = "slod", method = "hk", penalties = c(2.02, 2.62, 1.74) )
simspal <- calc.genoprob(simspal, step=0)
qtlslod <- makeqtl(simspal, chr = c(1, 4),
               pos = c(36.6, 27.8), what = "prob")


## ----profilelodimage, fig.height = 10, fig.width = 14--------------------
lodmat1.c <- getprofile(simspal, qtl =  qtlslod, pheno.cols =1:241,
                             formula = y~Q1 + Q2 , method = "hk",
							  verbose = F, tpy="comb")
plotprofile(lodmat1.c, mval = 8, col=heat.colors(100)[100:1], main="The Profile LOD image of data")


## ----effectplot----------------------------------------------------------

slodeff <- vector("list", nphe(simspal))

for(i in 1:nphe(simspal)) {
    slodeff[[i]] <- summary(fitqtl(simspal, phe=i, qtl=qtlslod,
                            method="hk", get.ests=TRUE, dropone=FALSE))$ests[,1]*c(1,2,2)
}

nam <- names(slodeff[[1]])
slodeff <- matrix(unlist(slodeff), byrow=TRUE, ncol=length(nam))
colnames(slodeff) <- nam

time <- (0:240)/30


## ----effplot, fig.height=5, fig.width = 14-------------------------------
par(mfrow=c(1,3))
plot(time, slodeff[,1], lwd=2, type="l",
     xlab="Time (hours)",
	      ylab="Tip angle (degrees)", col="red", ylim=c(-110,0))
		  mtext("baseline curve", side=3, line=0.5)

plot(time, slodeff[,2], lwd=2, ylim = c(-5,9), type="l",
     xlab="Time (hours)",
	      ylab="QTL effect", col="red")
		  abline(h=0)
		  mtext("chr 1, 37 cM", side=3, line=0.5)

plot(time, slodeff[,3], lwd=2, ylim = c(-5,9), type="l",
     xlab="Time (hours)",
	      ylab="QTL effect", col="red")
		  mtext("chr 4, 28 cM", side=3, line=0.5)
		  abline(h=0)


