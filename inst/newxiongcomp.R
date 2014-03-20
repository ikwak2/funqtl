## to get the simulation result described in Fig5 and 6

library(qtl)
library(fda)
source("fr.R")
source("logistic.R")
source("comparison_yap.R")
source("myfun20120901.R")
source("myfun20130116.R")


gen.data3 <- function(sample.size, cov.fcn, beta.coef,er){
  ##popu.size <- 10000  # population size

  ## simulate genotypes
  mp <- sim.map(100, n.mar=6, include.x=F, eq.spacing=T) # simulate map
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
#        D1 <- gen.data3(100, "autocorr", beta.coef,1 )


compsim3 <- function(iter, N, method, beta.coef,er ) {
    results = NULL;
    resultm = NULL;
    resultss = NULL;
    resultqf = NULL;
    for(rep in 1:iter) {
        D1 <- gen.data3(N, method, beta.coef,er)
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

permsim3 <- function(iter, N, method, beta.coef,er ) {
    results = NULL;
    resultm = NULL;
    resultss = NULL;
    resultqf = NULL;
    for(rep in 1:iter) {
        D1 <- gen.data3(N, method, beta.coef,er)
        o <- sample(N)

        pheno <- D1$pheno[o,];
        D1$pheno <-  pheno

        D1 <- calc.genoprob(D1, step=4)
        D1.out <- scanoneF(D1, pheno.cols=1:10, method="hk")

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

    Pos = reslod > lodcrit
    TP = ( abs(respos - 32) < window & reslod > lodcrit )
    FP = ( abs(respos - 32) >= window & reslod > lodcrit )
    FN = ( abs(respos - 32) < window & reslod <= lodcrit )
    TN = ( abs(respos - 32) >= window & reslod <= lodcrit )
    prtout <- list(length =  length(respos),
                   mean = mean(respos),
                   Pmean = mean(respos[reslod > lodcrit]),
                   sd = sd(respos),
                   rmse = sqrt(mean((respos-32)^2)),
                   TP = sum(TP), FP = sum(FP), FN = sum(FN), TN = sum(TN),
                   P = sum(Pos),
                   conf.int = binom.test(sum(Pos),10000)$conf.int
                   )
}

D1 <- gen.data3(200, "autocorr", beta.coef,2)
D2 <- gen.data3(200, "equicorr", beta.coef,2)
D3 <- gen.data3(200, "structured", beta.coef,2)

o <- sample(N)


# compxiongaut2.RData compxiongequi2.RData compxiongstr2.RData compxiong500aut2.RData compxiong500equi2.RData compxiong500str2.RData comp2xiong500aut2.RData comp2xiong500eq2.RData comp2xiong500st2.RData




out100aut1 <- compsim3(10000,100, "autocorr", beta.coef, 1)
out100eq1 <- compsim3(10000,100, "equicorr", beta.coef, 1)
out100st1 <- compsim3(10000,100, "structured", beta.coef, 1)

out100aut2 <- compsim3(10000,100, "autocorr", beta.coef, 2)
out100eq2 <- compsim3(10000,100, "equicorr", beta.coef, 2)
out100st2 <- compsim3(10000,100, "structured", beta.coef, 2)

out100aut3 <- compsim3(10000,100, "autocorr", beta.coef, 3)
out100eq3 <- compsim3(10000,100, "equicorr", beta.coef, 3)
out100st3 <- compsim3(10000,100, "structured", beta.coef, 3)

out100aut4 <- compsim3(10000,100, "autocorr", beta.coef, 6)
out100eq4 <- compsim3(10000,100, "equicorr", beta.coef, 6)
out100st4 <- compsim3(10000,100, "structured", beta.coef, .5)



out200aut1 <- compsim3(10000,200, "autocorr", beta.coef, 1)
out200eq1 <- compsim3(10000,200, "equicorr", beta.coef, 1)
out200st1 <- compsim3(10000,200, "structured", beta.coef, 1)

out200aut2 <- compsim3(10000,200, "autocorr", beta.coef, 2)
out200eq2 <- compsim3(10000,200, "equicorr", beta.coef, 2)
out200st2 <- compsim3(10000,200, "structured", beta.coef, 2)

out200aut3 <- compsim3(10000,200, "autocorr", beta.coef, 3)
out200eq3 <- compsim3(10000,200, "equicorr", beta.coef, 3)
out200st3 <- compsim3(10000,200, "structured", beta.coef, 3)

out200aut4 <- compsim3(10000,200, "autocorr", beta.coef, 6)
out200eq4 <- compsim3(10000,200, "equicorr", beta.coef, 6)
out200st4 <- compsim3(10000,200, "structured", beta.coef, .5)



out400aut1 <- compsim3(10000,400, "autocorr", beta.coef, 1)
out400eq1 <- compsim3(10000,400, "equicorr", beta.coef, 1)
out400st1 <- compsim3(10000,400, "structured", beta.coef, 1)

out400aut2 <- compsim3(10000,400, "autocorr", beta.coef, 2)
out400eq2 <- compsim3(10000,400, "equicorr", beta.coef, 2)
out400st2 <- compsim3(10000,400, "structured", beta.coef, 2)

out400aut3 <- compsim3(10000,400, "autocorr", beta.coef, 3)
out400eq3 <- compsim3(10000,400, "equicorr", beta.coef, 3)
out400st3 <- compsim3(10000,400, "structured", beta.coef, 3)

out400aut4 <- compsim3(10000,400, "autocorr", beta.coef, 6)
out400eq4 <- compsim3(10000,400, "equicorr", beta.coef, 6)
out400st4 <- compsim3(10000,400, "structured", beta.coef, .5)





pout100aut1 <- permsim3(10000,100, "autocorr", beta.coef, 1)
pout100eq1 <- permsim3(10000,100, "equicorr", beta.coef, 1)
pout100st1 <- permsim3(10000,100, "structured", beta.coef, 1)

pout100aut2 <- permsim3(10000,100, "autocorr", beta.coef, 2)
pout100eq2 <- permsim3(10000,100, "equicorr", beta.coef, 2)
pout100st2 <- permsim3(10000,100, "structured", beta.coef, 2)

pout100aut3 <- permsim3(10000,100, "autocorr", beta.coef, 3)
pout100eq3 <- permsim3(10000,100, "equicorr", beta.coef, 3)
pout100st3 <- permsim3(10000,100, "structured", beta.coef, 3)

pout100aut4 <- permsim3(10000,100, "autocorr", beta.coef, 6)
pout100eq4 <- permsim3(10000,100, "equicorr", beta.coef, 6)
pout100st4 <- permsim3(10000,100, "structured", beta.coef, .5)



pout200aut1 <- permsim3(10000,200, "autocorr", beta.coef, 1)
pout200eq1 <- permsim3(10000,200, "equicorr", beta.coef, 1)
pout200st1 <- permsim3(10000,200, "structured", beta.coef, 1)

pout200aut2 <- permsim3(10000,200, "autocorr", beta.coef, 2)
pout200eq2 <- permsim3(10000,200, "equicorr", beta.coef, 2)
pout200st2 <- permsim3(10000,200, "structured", beta.coef, 2)

pout200aut3 <- permsim3(10000,200, "autocorr", beta.coef, 3)
pout200eq3 <- permsim3(10000,200, "equicorr", beta.coef, 3)
pout200st3 <- permsim3(10000,200, "structured", beta.coef, 3)

pout200aut4 <- permsim3(10000,200, "autocorr", beta.coef, 6)
pout200eq4 <- permsim3(10000,200, "equicorr", beta.coef, 6)
pout200st4 <- permsim3(10000,200, "structured", beta.coef, .5)



pout400aut1 <- permsim3(10000,400, "autocorr", beta.coef, 1)
pout400eq1 <- permsim3(10000,400, "equicorr", beta.coef, 1)
pout400st1 <- permsim3(10000,400, "structured", beta.coef, 1)

pout400aut2 <- permsim3(10000,400, "autocorr", beta.coef, 2)
pout400eq2 <- permsim3(10000,400, "equicorr", beta.coef, 2)
pout400st2 <- permsim3(10000,400, "structured", beta.coef, 2)

pout400aut3 <- permsim3(10000,400, "autocorr", beta.coef, 3)
pout400eq3 <- permsim3(10000,400, "equicorr", beta.coef, 3)
pout400st3 <- permsim3(10000,400, "structured", beta.coef, 3)

pout400aut4 <- permsim3(10000,400, "autocorr", beta.coef, 6)
pout400eq4 <- permsim3(10000,400, "equicorr", beta.coef, 6)
pout400st4 <- permsim3(10000,400, "structured", beta.coef, .5)





load("outs1-1.RData")
load("outs1-2.RData")
load("outs2-1.RData")
load("outs2-2-1.RData")
load("outs2-2-2.RData")
load("outs3-1.RData")
load("outs3-2-1.RData")
load("outs3-2-2.RData")
load("outs4-1.RData")
load("outs4-2-1.RData")
load("outs4-2-2.RData")
load("outs1add200.RData")
load("outs2add200.RData")
load("outs3add200.RData")
load("outs4add200.RData")

out.100.au.1 <- out100aut1
out.100.eq.1 <- out100eq1
out.100.st.1 <- out100st1
out.400.au.1 <- out400aut1
out.400.eq.1 <- out400eq1
out.400.st.1 <- out400st1

out.100.au.2 <- out100aut2
out.100.eq.2 <- out100eq2
out.100.st.2 <- out100st2
out.400.au.2 <- out400aut2
out.400.eq.2 <- out400eq2
out.400.st.2 <- out400st2

out.100.au.3 <- out100aut3
out.100.eq.3 <- out100eq3
out.100.st.3 <- out100st3
out.400.au.3 <- out400aut3
out.400.eq.3 <- out400eq3
out.400.st.3 <- out400st3

out.100.au.4 <- out100aut4
out.100.eq.4 <- out100eq4
out.100.st.4 <- out100st4
out.400.au.4 <- out400aut4
out.400.eq.4 <- out400eq4
out.400.st.4 <- out400st4

out.200.au.1 <- out200aut1
out.200.eq.1 <- out200eq1
out.200.st.1 <- out200st1
out.200.au.2 <- out200aut2
out.200.eq.2 <- out200eq2
out.200.st.2 <- out200st2
out.200.au.3 <- out200aut3
out.200.eq.3 <- out200eq3
out.200.st.3 <- out200st3
out.200.au.4 <- out200aut4
out.200.eq.4 <- out200eq4
out.200.st.4 <- out200st4


#str(out.100.au.1)

load("pouts1-1-1.RData")
load("pouts1-1-2.RData")
load("pouts1-1-3.RData")
load("pouts1-2-1.RData")
load("pouts1-2-2.RData")
load("pouts1-2-3.RData")
load("pouts2-1-1.RData")
load("pouts2-1-2.RData")
load("pouts2-1-3.RData")
load("pouts2-2-1.RData")
load("pouts2-2-2.RData")
load("pouts2-2-3.RData")
load("pouts3-1-1.RData")
load("pouts3-1-2.RData")
load("pouts3-1-3.RData")
load("pouts3-2-1.RData")
load("pouts3-2-2.RData")
load("pouts3-2-3.RData")
load("pouts4-1-1.RData")
load("pouts4-1-2.RData")
load("pouts4-1-3.RData")
load("pouts4-2-1.RData")
load("pouts4-2-2.RData")
load("pouts4-2-3.RData")
load("poutsn1-1.RData")
load("poutsn1-2.RData")
load("poutsn1-3.RData")
load("poutsn2-1.RData")
load("poutsn2-2.RData")
load("poutsn2-3.RData")
load("poutsn3-1.RData")
load("poutsn3-2.RData")
load("poutsn3-3.RData")
load("poutsn4-1.RData")
load("poutsn4-2.RData")
load("poutsn4-3.RData")





(slod.aut.100.1 <- quantile(as.vector(pout100aut1$slod[,3]), 0.95))
(mlod.aut.100.1 <- quantile(as.vector(pout100aut1$mlod[,3]), 0.95))
(ss.aut.100.1 <- quantile(as.vector(pout100aut1$ss[,3]), 0.95))
(qf.aut.100.1 <- quantile(as.vector(pout100aut1$qf[,3]), 0.95))

(slod.aut.100.2 <- quantile(as.vector(pout100aut2$slod[,3]), 0.95))
(mlod.aut.100.2 <- quantile(as.vector(pout100aut2$mlod[,3]), 0.95))
(ss.aut.100.2 <- quantile(as.vector(pout100aut2$ss[,3]), 0.95))
(qf.aut.100.2 <- quantile(as.vector(pout100aut2$qf[,3]), 0.95))

(slod.aut.100.3 <- quantile(as.vector(pout100aut3$slod[,3]), 0.95))
(mlod.aut.100.3 <- quantile(as.vector(pout100aut3$mlod[,3]), 0.95))
(ss.aut.100.3 <- quantile(as.vector(pout100aut3$ss[,3]), 0.95))
(qf.aut.100.3 <- quantile(as.vector(pout100aut3$qf[,3]), 0.95))

(slod.aut.100.4 <- quantile(as.vector(pout100aut4$slod[,3]), 0.95))
(mlod.aut.100.4 <- quantile(as.vector(pout100aut4$mlod[,3]), 0.95))
(ss.aut.100.4 <- quantile(as.vector(pout100aut4$ss[,3]), 0.95))
(qf.aut.100.4 <- quantile(as.vector(pout100aut4$qf[,3]), 0.95))

(slod.aut.400.1 <- quantile(as.vector(pout400aut1$slod[,3]), 0.95))
(mlod.aut.400.1 <- quantile(as.vector(pout400aut1$mlod[,3]), 0.95))
(ss.aut.400.1 <- quantile(as.vector(pout400aut1$ss[,3]), 0.95))
(qf.aut.400.1 <- quantile(as.vector(pout400aut1$qf[,3]), 0.95))

(slod.aut.400.2 <- quantile(as.vector(pout400aut2$slod[,3]), 0.95))
(mlod.aut.400.2 <- quantile(as.vector(pout400aut2$mlod[,3]), 0.95))
(ss.aut.400.2 <- quantile(as.vector(pout400aut2$ss[,3]), 0.95))
(qf.aut.400.2 <- quantile(as.vector(pout400aut2$qf[,3]), 0.95))

(slod.aut.400.3 <- quantile(as.vector(pout400aut3$slod[,3]), 0.95))
(mlod.aut.400.3 <- quantile(as.vector(pout400aut3$mlod[,3]), 0.95))
(ss.aut.400.3 <- quantile(as.vector(pout400aut3$ss[,3]), 0.95))
(qf.aut.400.3 <- quantile(as.vector(pout400aut3$qf[,3]), 0.95))

(slod.aut.400.4 <- quantile(as.vector(pout400aut4$slod[,3]), 0.95))
(mlod.aut.400.4 <- quantile(as.vector(pout400aut4$mlod[,3]), 0.95))
(ss.aut.400.4 <- quantile(as.vector(pout400aut4$ss[,3]), 0.95))
(qf.aut.400.4 <- quantile(as.vector(pout400aut4$qf[,3]), 0.95))



(slod.eq.100.1 <- quantile(as.vector(pout100eq1$slod[,3]), 0.95))
(mlod.eq.100.1 <- quantile(as.vector(pout100eq1$mlod[,3]), 0.95))
(ss.eq.100.1 <- quantile(as.vector(pout100eq1$ss[,3]), 0.95))
(qf.eq.100.1 <- quantile(as.vector(pout100eq1$qf[,3]), 0.95))

(slod.eq.100.2 <- quantile(as.vector(pout100eq2$slod[,3]), 0.95))
(mlod.eq.100.2 <- quantile(as.vector(pout100eq2$mlod[,3]), 0.95))
(ss.eq.100.2 <- quantile(as.vector(pout100eq2$ss[,3]), 0.95))
(qf.eq.100.2 <- quantile(as.vector(pout100eq2$qf[,3]), 0.95))

(slod.eq.100.3 <- quantile(as.vector(pout100eq3$slod[,3]), 0.95))
(mlod.eq.100.3 <- quantile(as.vector(pout100eq3$mlod[,3]), 0.95))
(ss.eq.100.3 <- quantile(as.vector(pout100eq3$ss[,3]), 0.95))
(qf.eq.100.3 <- quantile(as.vector(pout100eq3$qf[,3]), 0.95))

(slod.eq.100.4 <- quantile(as.vector(pout100eq4$slod[,3]), 0.95))
(mlod.eq.100.4 <- quantile(as.vector(pout100eq4$mlod[,3]), 0.95))
(ss.eq.100.4 <- quantile(as.vector(pout100eq4$ss[,3]), 0.95))
(qf.eq.100.4 <- quantile(as.vector(pout100eq4$qf[,3]), 0.95))

(slod.eq.400.1 <- quantile(as.vector(pout400eq1$slod[,3]), 0.95))
(mlod.eq.400.1 <- quantile(as.vector(pout400eq1$mlod[,3]), 0.95))
(ss.eq.400.1 <- quantile(as.vector(pout400eq1$ss[,3]), 0.95))
(qf.eq.400.1 <- quantile(as.vector(pout400eq1$qf[,3]), 0.95))

(slod.eq.400.2 <- quantile(as.vector(pout400eq2$slod[,3]), 0.95))
(mlod.eq.400.2 <- quantile(as.vector(pout400eq2$mlod[,3]), 0.95))
(ss.eq.400.2 <- quantile(as.vector(pout400eq2$ss[,3]), 0.95))
(qf.eq.400.2 <- quantile(as.vector(pout400eq2$qf[,3]), 0.95))

(slod.eq.400.3 <- quantile(as.vector(pout400eq3$slod[,3]), 0.95))
(mlod.eq.400.3 <- quantile(as.vector(pout400eq3$mlod[,3]), 0.95))
(ss.eq.400.3 <- quantile(as.vector(pout400eq3$ss[,3]), 0.95))
(qf.eq.400.3 <- quantile(as.vector(pout400eq3$qf[,3]), 0.95))

(slod.eq.400.4 <- quantile(as.vector(pout400eq4$slod[,3]), 0.95))
(mlod.eq.400.4 <- quantile(as.vector(pout400eq4$mlod[,3]), 0.95))
(ss.eq.400.4 <- quantile(as.vector(pout400eq4$ss[,3]), 0.95))
(qf.eq.400.4 <- quantile(as.vector(pout400eq4$qf[,3]), 0.95))



(slod.st.100.1 <- quantile(as.vector(pout100st1$slod[,3]), 0.95))
(mlod.st.100.1 <- quantile(as.vector(pout100st1$mlod[,3]), 0.95))
(ss.st.100.1 <- quantile(as.vector(pout100st1$ss[,3]), 0.95))
(qf.st.100.1 <- quantile(as.vector(pout100st1$qf[,3]), 0.95))

(slod.st.100.2 <- quantile(as.vector(pout100st2$slod[,3]), 0.95))
(mlod.st.100.2 <- quantile(as.vector(pout100st2$mlod[,3]), 0.95))
(ss.st.100.2 <- quantile(as.vector(pout100st2$ss[,3]), 0.95))
(qf.st.100.2 <- quantile(as.vector(pout100st2$qf[,3]), 0.95))

(slod.st.100.3 <- quantile(as.vector(pout100st3$slod[,3]), 0.95))
(mlod.st.100.3 <- quantile(as.vector(pout100st3$mlod[,3]), 0.95))
(ss.st.100.3 <- quantile(as.vector(pout100st3$ss[,3]), 0.95))
(qf.st.100.3 <- quantile(as.vector(pout100st3$qf[,3]), 0.95))

(slod.st.100.4 <- quantile(as.vector(pout100st4$slod[,3]), 0.95))
(mlod.st.100.4 <- quantile(as.vector(pout100st4$mlod[,3]), 0.95))
(ss.st.100.4 <- quantile(as.vector(pout100st4$ss[,3]), 0.95))
(qf.st.100.4 <- quantile(as.vector(pout100st4$qf[,3]), 0.95))

(slod.st.400.1 <- quantile(as.vector(pout400st1$slod[,3]), 0.95))
(mlod.st.400.1 <- quantile(as.vector(pout400st1$mlod[,3]), 0.95))
(ss.st.400.1 <- quantile(as.vector(pout400st1$ss[,3]), 0.95))
(qf.st.400.1 <- quantile(as.vector(pout400st1$qf[,3]), 0.95))

(slod.st.400.2 <- quantile(as.vector(pout400st2$slod[,3]), 0.95))
(mlod.st.400.2 <- quantile(as.vector(pout400st2$mlod[,3]), 0.95))
(ss.st.400.2 <- quantile(as.vector(pout400st2$ss[,3]), 0.95))
(qf.st.400.2 <- quantile(as.vector(pout400st2$qf[,3]), 0.95))

(slod.st.400.3 <- quantile(as.vector(pout400st3$slod[,3]), 0.95))
(mlod.st.400.3 <- quantile(as.vector(pout400st3$mlod[,3]), 0.95))
(ss.st.400.3 <- quantile(as.vector(pout400st3$ss[,3]), 0.95))
(qf.st.400.3 <- quantile(as.vector(pout400st3$qf[,3]), 0.95))

(slod.st.400.4 <- quantile(as.vector(pout400st4$slod[,3]), 0.95))
(mlod.st.400.4 <- quantile(as.vector(pout400st4$mlod[,3]), 0.95))
(ss.st.400.4 <- quantile(as.vector(pout400st4$ss[,3]), 0.95))
(qf.st.400.4 <- quantile(as.vector(pout400st4$qf[,3]), 0.95))








(slod.aut.200.1 <- quantile(as.vector(pout200aut1$slod[,3]), 0.95))
(mlod.aut.200.1 <- quantile(as.vector(pout200aut1$mlod[,3]), 0.95))
(ss.aut.200.1 <- quantile(as.vector(pout200aut1$ss[,3]), 0.95))
(qf.aut.200.1 <- quantile(as.vector(pout200aut1$qf[,3]), 0.95))

(slod.aut.200.2 <- quantile(as.vector(pout200aut2$slod[,3]), 0.95))
(mlod.aut.200.2 <- quantile(as.vector(pout200aut2$mlod[,3]), 0.95))
(ss.aut.200.2 <- quantile(as.vector(pout200aut2$ss[,3]), 0.95))
(qf.aut.200.2 <- quantile(as.vector(pout200aut2$qf[,3]), 0.95))

(slod.aut.200.3 <- quantile(as.vector(pout200aut3$slod[,3]), 0.95))
(mlod.aut.200.3 <- quantile(as.vector(pout200aut3$mlod[,3]), 0.95))
(ss.aut.200.3 <- quantile(as.vector(pout200aut3$ss[,3]), 0.95))
(qf.aut.200.3 <- quantile(as.vector(pout200aut3$qf[,3]), 0.95))

(slod.aut.200.4 <- quantile(as.vector(pout200aut4$slod[,3]), 0.95))
(mlod.aut.200.4 <- quantile(as.vector(pout200aut4$mlod[,3]), 0.95))
(ss.aut.200.4 <- quantile(as.vector(pout200aut4$ss[,3]), 0.95))
(qf.aut.200.4 <- quantile(as.vector(pout200aut4$qf[,3]), 0.95))

(slod.eq.200.1 <- quantile(as.vector(pout200eq1$slod[,3]), 0.95))
(mlod.eq.200.1 <- quantile(as.vector(pout200eq1$mlod[,3]), 0.95))
(ss.eq.200.1 <- quantile(as.vector(pout200eq1$ss[,3]), 0.95))
(qf.eq.200.1 <- quantile(as.vector(pout200eq1$qf[,3]), 0.95))

(slod.eq.200.2 <- quantile(as.vector(pout200eq2$slod[,3]), 0.95))
(mlod.eq.200.2 <- quantile(as.vector(pout200eq2$mlod[,3]), 0.95))
(ss.eq.200.2 <- quantile(as.vector(pout200eq2$ss[,3]), 0.95))
(qf.eq.200.2 <- quantile(as.vector(pout200eq2$qf[,3]), 0.95))

(slod.eq.200.3 <- quantile(as.vector(pout200eq3$slod[,3]), 0.95))
(mlod.eq.200.3 <- quantile(as.vector(pout200eq3$mlod[,3]), 0.95))
(ss.eq.200.3 <- quantile(as.vector(pout200eq3$ss[,3]), 0.95))
(qf.eq.200.3 <- quantile(as.vector(pout200eq3$qf[,3]), 0.95))

(slod.eq.200.4 <- quantile(as.vector(pout200eq4$slod[,3]), 0.95))
(mlod.eq.200.4 <- quantile(as.vector(pout200eq4$mlod[,3]), 0.95))
(ss.eq.200.4 <- quantile(as.vector(pout200eq4$ss[,3]), 0.95))
(qf.eq.200.4 <- quantile(as.vector(pout200eq4$qf[,3]), 0.95))

(slod.st.200.1 <- quantile(as.vector(pout200st1$slod[,3]), 0.95))
(mlod.st.200.1 <- quantile(as.vector(pout200st1$mlod[,3]), 0.95))
(ss.st.200.1 <- quantile(as.vector(pout200st1$ss[,3]), 0.95))
(qf.st.200.1 <- quantile(as.vector(pout200st1$qf[,3]), 0.95))

(slod.st.200.2 <- quantile(as.vector(pout200st2$slod[,3]), 0.95))
(mlod.st.200.2 <- quantile(as.vector(pout200st2$mlod[,3]), 0.95))
(ss.st.200.2 <- quantile(as.vector(pout200st2$ss[,3]), 0.95))
(qf.st.200.2 <- quantile(as.vector(pout200st2$qf[,3]), 0.95))

(slod.st.200.3 <- quantile(as.vector(pout200st3$slod[,3]), 0.95))
(mlod.st.200.3 <- quantile(as.vector(pout200st3$mlod[,3]), 0.95))
(ss.st.200.3 <- quantile(as.vector(pout200st3$ss[,3]), 0.95))
(qf.st.200.3 <- quantile(as.vector(pout200st3$qf[,3]), 0.95))

(slod.st.200.4 <- quantile(as.vector(pout200st4$slod[,3]), 0.95))
(mlod.st.200.4 <- quantile(as.vector(pout200st4$mlod[,3]), 0.95))
(ss.st.200.4 <- quantile(as.vector(pout200st4$ss[,3]), 0.95))
(qf.st.200.4 <- quantile(as.vector(pout200st4$qf[,3]), 0.95))











# out100aut1

(pout.100.au.1.s <- prtout(out.100.au.1$slod, slod.aut.100.1, 7))
(pout.100.au.1.m <- prtout(out.100.au.1$mlod, mlod.aut.100.1, 7))
(pout.100.au.1.ss <- prtout(out.100.au.1$ss, ss.aut.100.1, 7))
(pout.100.au.1.qf <- prtout(out.100.au.1$qf, qf.aut.100.1, 7))

(pout.100.au.2.s <- prtout(out.100.au.2$slod, slod.aut.100.2, 7))
(pout.100.au.2.m <- prtout(out.100.au.2$mlod, mlod.aut.100.2, 7))
(pout.100.au.2.ss <- prtout(out.100.au.2$ss, ss.aut.100.2, 7))
(pout.100.au.2.qf <- prtout(out.100.au.2$qf, qf.aut.100.2, 7))

(pout.100.au.3.s <- prtout(out.100.au.3$slod, slod.aut.100.3, 7))
(pout.100.au.3.m <- prtout(out.100.au.3$mlod, mlod.aut.100.3, 7))
(pout.100.au.3.ss <- prtout(out.100.au.3$ss, ss.aut.100.3, 7))
(pout.100.au.3.qf <- prtout(out.100.au.3$qf, qf.aut.100.3, 7))

(pout.100.au.4.s <- prtout(out.100.au.4$slod, slod.aut.100.4, 7))
(pout.100.au.4.m <- prtout(out.100.au.4$mlod, mlod.aut.100.4, 7))
(pout.100.au.4.ss <- prtout(out.100.au.4$ss, ss.aut.100.4, 7))
(pout.100.au.4.qf <- prtout(out.100.au.4$qf, qf.aut.100.4, 7))




(pout.100.eq.1.s <- prtout(out.100.eq.1$slod, slod.eq.100.1, 7))
(pout.100.eq.1.m <- prtout(out.100.eq.1$mlod, mlod.eq.100.1, 7))
(pout.100.eq.1.ss <- prtout(out.100.eq.1$ss, ss.eq.100.1, 7))
(pout.100.eq.1.qf <- prtout(out.100.eq.1$qf, qf.eq.100.1, 7))

(pout.100.eq.2.s <- prtout(out.100.eq.2$slod, slod.eq.100.2, 7))
(pout.100.eq.2.m <- prtout(out.100.eq.2$mlod, mlod.eq.100.2, 7))
(pout.100.eq.2.ss <- prtout(out.100.eq.2$ss, ss.eq.100.2, 7))
(pout.100.eq.2.qf <- prtout(out.100.eq.2$qf, qf.eq.100.2, 7))

(pout.100.eq.3.s <- prtout(out.100.eq.3$slod, slod.eq.100.3, 7))
(pout.100.eq.3.m <- prtout(out.100.eq.3$mlod, mlod.eq.100.3, 7))
(pout.100.eq.3.ss <- prtout(out.100.eq.3$ss, ss.eq.100.3, 7))
(pout.100.eq.3.qf <- prtout(out.100.eq.3$qf, qf.eq.100.3, 7))

(pout.100.eq.4.s <- prtout(out.100.eq.4$slod, slod.eq.100.4, 7))
(pout.100.eq.4.m <- prtout(out.100.eq.4$mlod, mlod.eq.100.4, 7))
(pout.100.eq.4.ss <- prtout(out.100.eq.4$ss, ss.eq.100.4, 7))
(pout.100.eq.4.qf <- prtout(out.100.eq.4$qf, qf.eq.100.4, 7))








(pout.100.st.1.s <- prtout(out.100.st.1$slod, slod.st.100.1, 7))
(pout.100.st.1.m <- prtout(out.100.st.1$mlod, mlod.st.100.1, 7))
(pout.100.st.1.ss <- prtout(out.100.st.1$ss, ss.st.100.1, 7))
(pout.100.st.1.qf <- prtout(out.100.st.1$qf, qf.st.100.1, 7))

(pout.100.st.2.s <- prtout(out.100.st.2$slod, slod.st.100.2, 7))
(pout.100.st.2.m <- prtout(out.100.st.2$mlod, mlod.st.100.2, 7))
(pout.100.st.2.ss <- prtout(out.100.st.2$ss, ss.st.100.2, 7))
(pout.100.st.2.qf <- prtout(out.100.st.2$qf, qf.st.100.2, 7))

(pout.100.st.3.s <- prtout(out.100.st.3$slod, slod.st.100.3, 7))
(pout.100.st.3.m <- prtout(out.100.st.3$mlod, mlod.st.100.3, 7))
(pout.100.st.3.ss <- prtout(out.100.st.3$ss, ss.st.100.3, 7))
(pout.100.st.3.qf <- prtout(out.100.st.3$qf, qf.st.100.3, 7))

(pout.100.st.4.s <- prtout(out.100.st.4$slod, slod.st.100.4, 7))
(pout.100.st.4.m <- prtout(out.100.st.4$mlod, mlod.st.100.4, 7))
(pout.100.st.4.ss <- prtout(out.100.st.4$ss, ss.st.100.4, 7))
(pout.100.st.4.qf <- prtout(out.100.st.4$qf, qf.st.100.4, 7))




(pout.400.au.1.s <- prtout(out.400.au.1$slod, slod.aut.400.1, 7))
(pout.400.au.1.m <- prtout(out.400.au.1$mlod, mlod.aut.400.1, 7))
(pout.400.au.1.ss <- prtout(out.400.au.1$ss, ss.aut.400.1, 7))
(pout.400.au.1.qf <- prtout(out.400.au.1$qf, qf.aut.400.1, 7))

(pout.400.au.2.s <- prtout(out.400.au.2$slod, slod.aut.400.2, 7))
(pout.400.au.2.m <- prtout(out.400.au.2$mlod, mlod.aut.400.2, 7))
(pout.400.au.2.ss <- prtout(out.400.au.2$ss, ss.aut.400.2, 7))
(pout.400.au.2.qf <- prtout(out.400.au.2$qf, qf.aut.400.2, 7))

(pout.400.au.3.s <- prtout(out.400.au.3$slod, slod.aut.400.3, 7))
(pout.400.au.3.m <- prtout(out.400.au.3$mlod, mlod.aut.400.3, 7))
(pout.400.au.3.ss <- prtout(out.400.au.3$ss, ss.aut.400.3, 7))
(pout.400.au.3.qf <- prtout(out.400.au.3$qf, qf.aut.400.3, 7))

(pout.400.au.4.s <- prtout(out.400.au.4$slod, slod.aut.400.4, 7))
(pout.400.au.4.m <- prtout(out.400.au.4$mlod, mlod.aut.400.4, 7))
(pout.400.au.4.ss <- prtout(out.400.au.4$ss, ss.aut.400.4, 7))
(pout.400.au.4.qf <- prtout(out.400.au.4$qf, qf.aut.400.4, 7))





(pout.400.eq.1.s <- prtout(out.400.eq.1$slod, slod.eq.400.1, 7))
(pout.400.eq.1.m <- prtout(out.400.eq.1$mlod, mlod.eq.400.1, 7))
(pout.400.eq.1.ss <- prtout(out.400.eq.1$ss, ss.eq.400.1, 7))
(pout.400.eq.1.qf <- prtout(out.400.eq.1$qf, qf.eq.400.1, 7))

(pout.400.eq.2.s <- prtout(out.400.eq.2$slod, slod.eq.400.2, 7))
(pout.400.eq.2.m <- prtout(out.400.eq.2$mlod, mlod.eq.400.2, 7))
(pout.400.eq.2.ss <- prtout(out.400.eq.2$ss, ss.eq.400.2, 7))
(pout.400.eq.2.qf <- prtout(out.400.eq.2$qf, qf.eq.400.2, 7))

(pout.400.eq.3.s <- prtout(out.400.eq.3$slod, slod.eq.400.3, 7))
(pout.400.eq.3.m <- prtout(out.400.eq.3$mlod, mlod.eq.400.3, 7))
(pout.400.eq.3.ss <- prtout(out.400.eq.3$ss, ss.eq.400.3, 7))
(pout.400.eq.3.qf <- prtout(out.400.eq.3$qf, qf.eq.400.3, 7))

(pout.400.eq.4.s <- prtout(out.400.eq.4$slod, slod.eq.400.4, 7))
(pout.400.eq.4.m <- prtout(out.400.eq.4$mlod, mlod.eq.400.4, 7))
(pout.400.eq.4.ss <- prtout(out.400.eq.4$ss, ss.eq.400.4, 7))
(pout.400.eq.4.qf <- prtout(out.400.eq.4$qf, qf.eq.400.4, 7))



(pout.400.st.1.s <- prtout(out.400.st.1$slod, slod.st.400.1, 7))
(pout.400.st.1.m <- prtout(out.400.st.1$mlod, mlod.st.400.1, 7))
(pout.400.st.1.ss <- prtout(out.400.st.1$ss, ss.st.400.1, 7))
(pout.400.st.1.qf <- prtout(out.400.st.1$qf, qf.st.400.1, 7))

(pout.400.st.2.s <- prtout(out.400.st.2$slod, slod.st.400.2, 7))
(pout.400.st.2.m <- prtout(out.400.st.2$mlod, mlod.st.400.2, 7))
(pout.400.st.2.ss <- prtout(out.400.st.2$ss, ss.st.400.2, 7))
(pout.400.st.2.qf <- prtout(out.400.st.2$qf, qf.st.400.2, 7))

(pout.400.st.3.s <- prtout(out.400.st.3$slod, slod.st.400.3, 7))
(pout.400.st.3.m <- prtout(out.400.st.3$mlod, mlod.st.400.3, 7))
(pout.400.st.3.ss <- prtout(out.400.st.3$ss, ss.st.400.3, 7))
(pout.400.st.3.qf <- prtout(out.400.st.3$qf, qf.st.400.3, 7))

(pout.400.st.4.s <- prtout(out.400.st.4$slod, slod.st.400.4, 7))
(pout.400.st.4.m <- prtout(out.400.st.4$mlod, mlod.st.400.4, 7))
(pout.400.st.4.ss <- prtout(out.400.st.4$ss, ss.st.400.4, 7))
(pout.400.st.4.qf <- prtout(out.400.st.4$qf, qf.st.400.4, 7))





(pout.200.au.1.s <- prtout(out.200.au.1$slod, slod.aut.200.1, 7))
(pout.200.au.1.m <- prtout(out.200.au.1$mlod, mlod.aut.200.1, 7))
(pout.200.au.1.ss <- prtout(out.200.au.1$ss, ss.aut.200.1, 7))
(pout.200.au.1.qf <- prtout(out.200.au.1$qf, qf.aut.200.1, 7))

(pout.200.au.2.s <- prtout(out.200.au.2$slod, slod.aut.200.2, 7))
(pout.200.au.2.m <- prtout(out.200.au.2$mlod, mlod.aut.200.2, 7))
(pout.200.au.2.ss <- prtout(out.200.au.2$ss, ss.aut.200.2, 7))
(pout.200.au.2.qf <- prtout(out.200.au.2$qf, qf.aut.200.2, 7))

(pout.200.au.3.s <- prtout(out.200.au.3$slod, slod.aut.200.3, 7))
(pout.200.au.3.m <- prtout(out.200.au.3$mlod, mlod.aut.200.3, 7))
(pout.200.au.3.ss <- prtout(out.200.au.3$ss, ss.aut.200.3, 7))
(pout.200.au.3.qf <- prtout(out.200.au.3$qf, qf.aut.200.3, 7))

(pout.200.au.4.s <- prtout(out.200.au.4$slod, slod.aut.200.4, 7))
(pout.200.au.4.m <- prtout(out.200.au.4$mlod, mlod.aut.200.4, 7))
(pout.200.au.4.ss <- prtout(out.200.au.4$ss, ss.aut.200.4, 7))
(pout.200.au.4.qf <- prtout(out.200.au.4$qf, qf.aut.200.4, 7))




(pout.200.eq.1.s <- prtout(out.200.eq.1$slod, slod.eq.200.1, 7))
(pout.200.eq.1.m <- prtout(out.200.eq.1$mlod, mlod.eq.200.1, 7))
(pout.200.eq.1.ss <- prtout(out.200.eq.1$ss, ss.eq.200.1, 7))
(pout.200.eq.1.qf <- prtout(out.200.eq.1$qf, qf.eq.200.1, 7))

(pout.200.eq.2.s <- prtout(out.200.eq.2$slod, slod.eq.200.2, 7))
(pout.200.eq.2.m <- prtout(out.200.eq.2$mlod, mlod.eq.200.2, 7))
(pout.200.eq.2.ss <- prtout(out.200.eq.2$ss, ss.eq.200.2, 7))
(pout.200.eq.2.qf <- prtout(out.200.eq.2$qf, qf.eq.200.2, 7))

(pout.200.eq.3.s <- prtout(out.200.eq.3$slod, slod.eq.200.3, 7))
(pout.200.eq.3.m <- prtout(out.200.eq.3$mlod, mlod.eq.200.3, 7))
(pout.200.eq.3.ss <- prtout(out.200.eq.3$ss, ss.eq.200.3, 7))
(pout.200.eq.3.qf <- prtout(out.200.eq.3$qf, qf.eq.200.3, 7))

(pout.200.eq.4.s <- prtout(out.200.eq.4$slod, slod.eq.200.4, 7))
(pout.200.eq.4.m <- prtout(out.200.eq.4$mlod, mlod.eq.200.4, 7))
(pout.200.eq.4.ss <- prtout(out.200.eq.4$ss, ss.eq.200.4, 7))
(pout.200.eq.4.qf <- prtout(out.200.eq.4$qf, qf.eq.200.4, 7))








(pout.200.st.1.s <- prtout(out.200.st.1$slod, slod.st.200.1, 7))
(pout.200.st.1.m <- prtout(out.200.st.1$mlod, mlod.st.200.1, 7))
(pout.200.st.1.ss <- prtout(out.200.st.1$ss, ss.st.200.1, 7))
(pout.200.st.1.qf <- prtout(out.200.st.1$qf, qf.st.200.1, 7))

(pout.200.st.2.s <- prtout(out.200.st.2$slod, slod.st.200.2, 7))
(pout.200.st.2.m <- prtout(out.200.st.2$mlod, mlod.st.200.2, 7))
(pout.200.st.2.ss <- prtout(out.200.st.2$ss, ss.st.200.2, 7))
(pout.200.st.2.qf <- prtout(out.200.st.2$qf, qf.st.200.2, 7))

(pout.200.st.3.s <- prtout(out.200.st.3$slod, slod.st.200.3, 7))
(pout.200.st.3.m <- prtout(out.200.st.3$mlod, mlod.st.200.3, 7))
(pout.200.st.3.ss <- prtout(out.200.st.3$ss, ss.st.200.3, 7))
(pout.200.st.3.qf <- prtout(out.200.st.3$qf, qf.st.200.3, 7))

(pout.200.st.4.s <- prtout(out.200.st.4$slod, slod.st.200.4, 7))
(pout.200.st.4.m <- prtout(out.200.st.4$mlod, mlod.st.200.4, 7))
(pout.200.st.4.ss <- prtout(out.200.st.4$ss, ss.st.200.4, 7))
(pout.200.st.4.qf <- prtout(out.200.st.4$qf, qf.st.200.4, 7))






#### 100 - autocorr

aut100rmsea <- cbind(c(pout.100.au.1.s$rmse, pout.100.au.1.m$rmse, pout.100.au.1.ss$rmse, pout.100.au.1.qf$rmse),
                     c(pout.100.au.2.s$rmse, pout.100.au.2.m$rmse, pout.100.au.2.ss$rmse, pout.100.au.2.qf$rmse),
                     c(pout.100.au.3.s$rmse, pout.100.au.3.m$rmse, pout.100.au.3.ss$rmse, pout.100.au.3.qf$rmse),
                     c(pout.100.au.4.s$rmse, pout.100.au.4.m$rmse, pout.100.au.4.ss$rmse, pout.100.au.4.qf$rmse)
                     )

aut100powersa <- cbind(c(pout.100.au.1.s$P, pout.100.au.1.m$P, pout.100.au.1.ss$P, pout.100.au.1.qf$P),
                       c(pout.100.au.2.s$P, pout.100.au.2.m$P, pout.100.au.2.ss$P, pout.100.au.2.qf$P),
                      c(pout.100.au.3.s$P, pout.100.au.3.m$P, pout.100.au.3.ss$P, pout.100.au.3.qf$P),
                      c(pout.100.au.4.s$P, pout.100.au.4.m$P, pout.100.au.4.ss$P, pout.100.au.4.qf$P)
                      ) / 10000

### 100 - equicorr

aut100rmsee <- cbind(c(pout.100.eq.1.s$rmse, pout.100.eq.1.m$rmse, pout.100.eq.1.ss$rmse, pout.100.eq.1.qf$rmse),
                    c(pout.100.eq.2.s$rmse, pout.100.eq.2.m$rmse, pout.100.eq.2.ss$rmse, pout.100.eq.2.qf$rmse),
                    c(pout.100.eq.3.s$rmse, pout.100.eq.3.m$rmse, pout.100.eq.3.ss$rmse, pout.100.eq.3.qf$rmse),
                    c(pout.100.eq.4.s$rmse, pout.100.eq.4.m$rmse, pout.100.eq.4.ss$rmse, pout.100.eq.4.qf$rmse)
                    )

aut100powerse <- cbind(c(pout.100.eq.1.s$P, pout.100.eq.1.m$P, pout.100.eq.1.ss$P, pout.100.eq.1.qf$P),
                      c(pout.100.eq.2.s$P, pout.100.eq.2.m$P, pout.100.eq.2.ss$P, pout.100.eq.2.qf$P),
                      c(pout.100.eq.3.s$P, pout.100.eq.3.m$P, pout.100.eq.3.ss$P, pout.100.eq.3.qf$P),
                      c(pout.100.eq.4.s$P, pout.100.eq.4.m$P, pout.100.eq.4.ss$P, pout.100.eq.4.qf$P)
                      ) / 10000


### 100 - structured

aut100rmses <- cbind(c(pout.100.st.4.s$rmse, pout.100.st.4.m$rmse, pout.100.st.4.ss$rmse, pout.100.st.4.qf$rmse),
                    c(pout.100.st.1.s$rmse, pout.100.st.1.m$rmse, pout.100.st.1.ss$rmse, pout.100.st.1.qf$rmse),
                    c(pout.100.st.2.s$rmse, pout.100.st.2.m$rmse, pout.100.st.2.ss$rmse, pout.100.st.2.qf$rmse),
                    c(pout.100.st.3.s$rmse, pout.100.st.3.m$rmse, pout.100.st.3.ss$rmse, pout.100.st.3.qf$rmse)
                    )

aut100powerss <- cbind(c(pout.100.st.4.s$P, pout.100.st.4.m$P, pout.100.st.4.ss$P, pout.100.st.4.qf$P),
                      c(pout.100.st.1.s$P, pout.100.st.1.m$P, pout.100.st.1.ss$P, pout.100.st.1.qf$P),
                      c(pout.100.st.2.s$P, pout.100.st.2.m$P, pout.100.st.2.ss$P, pout.100.st.2.qf$P),
                      c(pout.100.st.3.s$P, pout.100.st.3.m$P, pout.100.st.3.ss$P, pout.100.st.3.qf$P)
                      ) / 10000




#### 200 - autocorr

aut200rmsea <- cbind(c(pout.200.au.1.s$rmse, pout.200.au.1.m$rmse, pout.200.au.1.ss$rmse, pout.200.au.1.qf$rmse),
                     c(pout.200.au.2.s$rmse, pout.200.au.2.m$rmse, pout.200.au.2.ss$rmse, pout.200.au.2.qf$rmse),
                     c(pout.200.au.3.s$rmse, pout.200.au.3.m$rmse, pout.200.au.3.ss$rmse, pout.200.au.3.qf$rmse),
                     c(pout.200.au.4.s$rmse, pout.200.au.4.m$rmse, pout.200.au.4.ss$rmse, pout.200.au.4.qf$rmse)
                     )

aut200powersa <- cbind(c(pout.200.au.1.s$P, pout.200.au.1.m$P, pout.200.au.1.ss$P, pout.200.au.1.qf$P),
                       c(pout.200.au.2.s$P, pout.200.au.2.m$P, pout.200.au.2.ss$P, pout.200.au.2.qf$P),
                      c(pout.200.au.3.s$P, pout.200.au.3.m$P, pout.200.au.3.ss$P, pout.200.au.3.qf$P),
                      c(pout.200.au.4.s$P, pout.200.au.4.m$P, pout.200.au.4.ss$P, pout.200.au.4.qf$P)
                      ) / 10000

### 200 - equicorr

aut200rmsee <- cbind(c(pout.200.eq.1.s$rmse, pout.200.eq.1.m$rmse, pout.200.eq.1.ss$rmse, pout.200.eq.1.qf$rmse),
                    c(pout.200.eq.2.s$rmse, pout.200.eq.2.m$rmse, pout.200.eq.2.ss$rmse, pout.200.eq.2.qf$rmse),
                    c(pout.200.eq.3.s$rmse, pout.200.eq.3.m$rmse, pout.200.eq.3.ss$rmse, pout.200.eq.3.qf$rmse),
                    c(pout.200.eq.4.s$rmse, pout.200.eq.4.m$rmse, pout.200.eq.4.ss$rmse, pout.200.eq.4.qf$rmse)
                    )

aut200powerse <- cbind(c(pout.200.eq.1.s$P, pout.200.eq.1.m$P, pout.200.eq.1.ss$P, pout.200.eq.1.qf$P),
                      c(pout.200.eq.2.s$P, pout.200.eq.2.m$P, pout.200.eq.2.ss$P, pout.200.eq.2.qf$P),
                      c(pout.200.eq.3.s$P, pout.200.eq.3.m$P, pout.200.eq.3.ss$P, pout.200.eq.3.qf$P),
                      c(pout.200.eq.4.s$P, pout.200.eq.4.m$P, pout.200.eq.4.ss$P, pout.200.eq.4.qf$P)
                      ) / 10000


### 200 - structured

aut200rmses <- cbind(c(pout.200.st.4.s$rmse, pout.200.st.4.m$rmse, pout.200.st.4.ss$rmse, pout.200.st.4.qf$rmse),
                    c(pout.200.st.1.s$rmse, pout.200.st.1.m$rmse, pout.200.st.1.ss$rmse, pout.200.st.1.qf$rmse),
                    c(pout.200.st.2.s$rmse, pout.200.st.2.m$rmse, pout.200.st.2.ss$rmse, pout.200.st.2.qf$rmse),
                    c(pout.200.st.3.s$rmse, pout.200.st.3.m$rmse, pout.200.st.3.ss$rmse, pout.200.st.3.qf$rmse)
                    )

aut200powerss <- cbind(c(pout.200.st.4.s$P, pout.200.st.4.m$P, pout.200.st.4.ss$P, pout.200.st.4.qf$P),
                      c(pout.200.st.1.s$P, pout.200.st.1.m$P, pout.200.st.1.ss$P, pout.200.st.1.qf$P),
                      c(pout.200.st.2.s$P, pout.200.st.2.m$P, pout.200.st.2.ss$P, pout.200.st.2.qf$P),
                      c(pout.200.st.3.s$P, pout.200.st.3.m$P, pout.200.st.3.ss$P, pout.200.st.3.qf$P)
                      ) / 10000



### 400 - autocorr

aut400rmsea <- cbind(c(pout.400.au.1.s$rmse, pout.400.au.1.m$rmse, pout.400.au.1.ss$rmse, pout.400.au.1.qf$rmse),
                     c(pout.400.au.2.s$rmse, pout.400.au.2.m$rmse, pout.400.au.2.ss$rmse, pout.400.au.2.qf$rmse),
                     c(pout.400.au.3.s$rmse, pout.400.au.3.m$rmse, pout.400.au.3.ss$rmse, pout.400.au.3.qf$rmse),
                     c(pout.400.au.4.s$rmse, pout.400.au.4.m$rmse, pout.400.au.4.ss$rmse, pout.400.au.4.qf$rmse)
                    )

aut400powersa <- cbind(c(pout.400.au.1.s$P, pout.400.au.1.m$P, pout.400.au.1.ss$P, pout.400.au.1.qf$P),
                      c(pout.400.au.2.s$P, pout.400.au.2.m$P, pout.400.au.2.ss$P, pout.400.au.2.qf$P),
                      c(pout.400.au.3.s$P, pout.400.au.3.m$P, pout.400.au.3.ss$P, pout.400.au.3.qf$P),
                      c(pout.400.au.4.s$P, pout.400.au.4.m$P, pout.400.au.4.ss$P, pout.400.au.4.qf$P)
                      ) / 10000

### 400 - equicorr

aut400rmsee <- cbind(c(pout.400.eq.1.s$rmse, pout.400.eq.1.m$rmse, pout.400.eq.1.ss$rmse, pout.400.eq.1.qf$rmse),
                    c(pout.400.eq.2.s$rmse, pout.400.eq.2.m$rmse, pout.400.eq.2.ss$rmse, pout.400.eq.2.qf$rmse),
                    c(pout.400.eq.3.s$rmse, pout.400.eq.3.m$rmse, pout.400.eq.3.ss$rmse, pout.400.eq.3.qf$rmse),
                    c(pout.400.eq.4.s$rmse, pout.400.eq.4.m$rmse, pout.400.eq.4.ss$rmse, pout.400.eq.4.qf$rmse)
                    )

aut400powerse <- cbind(c(pout.400.eq.1.s$P, pout.400.eq.1.m$P, pout.400.eq.1.ss$P, pout.400.eq.1.qf$P),
                      c(pout.400.eq.2.s$P, pout.400.eq.2.m$P, pout.400.eq.2.ss$P, pout.400.eq.2.qf$P),
                      c(pout.400.eq.3.s$P, pout.400.eq.3.m$P, pout.400.eq.3.ss$P, pout.400.eq.3.qf$P),
                      c(pout.400.eq.4.s$P, pout.400.eq.4.m$P, pout.400.eq.4.ss$P, pout.400.eq.4.qf$P)
                      ) / 10000


### 400 - structured

aut400rmses <- cbind(c(pout.400.st.4.s$rmse, pout.400.st.4.m$rmse, pout.400.st.4.ss$rmse, pout.400.st.4.qf$rmse),
                    c(pout.400.st.1.s$rmse, pout.400.st.1.m$rmse, pout.400.st.1.ss$rmse, pout.400.st.1.qf$rmse),
                    c(pout.400.st.2.s$rmse, pout.400.st.2.m$rmse, pout.400.st.2.ss$rmse, pout.400.st.2.qf$rmse),
                    c(pout.400.st.3.s$rmse, pout.400.st.3.m$rmse, pout.400.st.3.ss$rmse, pout.400.st.3.qf$rmse)
                    )

aut400powerss <- cbind(c(pout.400.st.4.s$P, pout.400.st.4.m$P, pout.400.st.4.ss$P, pout.400.st.4.qf$P),
                      c(pout.400.st.1.s$P, pout.400.st.1.m$P, pout.400.st.1.ss$P, pout.400.st.1.qf$P),
                      c(pout.400.st.2.s$P, pout.400.st.2.m$P, pout.400.st.2.ss$P, pout.400.st.2.qf$P),
                      c(pout.400.st.3.s$P, pout.400.st.3.m$P, pout.400.st.3.ss$P, pout.400.st.3.qf$P)
                      ) / 10000



#####


t1 <- c(7.992522, 4.417868, 3.085633, 1.504872)
t2 <- c(8.391845, 4.371283, 3.111750, 1.569568)
t3 <- c(6.775272, 3.618401, 2.013146, 1.175318)

par(mfrow=c(3,3))

plot(t1, c(aut100powersa[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="aut100", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t1, c(aut100powersa[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t1, c(aut100powersa[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t1, c(aut100powersa[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")


plot(t2, c(aut100powerse[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="equ100", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t2, c(aut100powerse[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t2, c(aut100powerse[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t2, c(aut100powerse[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")


plot(t3, c(aut100powerss[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="st100", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t3, c(aut100powerss[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t3, c(aut100powerss[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t3, c(aut100powerss[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")



plot(t1, c(aut400powersa[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="aut400", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t1, c(aut400powersa[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t1, c(aut400powersa[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t1, c(aut400powersa[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")


plot(t2, c(aut400powerse[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="equ400", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t2, c(aut400powerse[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t2, c(aut400powerse[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t2, c(aut400powerse[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")


plot(t3, c(aut400powerss[1,]), xlim=c(0,10), ylim=c(0,1), type="l", col="red", main="st400", ylab="Power", xlab="Percent variance explained by QTL")
par(new=T);plot(t3, c(aut400powerss[2,]), xlim=c(0,10), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(t3, c(aut400powerss[3,]), xlim=c(0,10), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")
par(new=T);plot(t3, c(aut400powerss[4,]), xlim=c(0,10), ylim=c(0,1), type="l", col="brown", ylab="Power", xlab="")


par(mfrow=c(2,3))

plot(t1, c(aut100rmsea[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="aut100", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t1, c(aut100rmsea[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t1, c(aut100rmsea[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t1, c(aut100rmsea[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")

plot(t2, c(aut100rmsee[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="eq100", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t2, c(aut100rmsee[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t2, c(aut100rmsee[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t2, c(aut100rmsee[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")

plot(t3, c(aut100rmses[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="st100", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t3, c(aut100rmses[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t3, c(aut100rmses[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t3, c(aut100rmses[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")



plot(t1, c(aut400rmsea[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="aut400", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t1, c(aut400rmsea[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t1, c(aut400rmsea[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t1, c(aut400rmsea[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")

plot(t2, c(aut400rmsee[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="eq400", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t2, c(aut400rmsee[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t2, c(aut400rmsee[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t2, c(aut400rmsee[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")

plot(t3, c(aut400rmses[1,]), xlim=c(0,10), ylim=c(0,40), type="l", col="red", main="st400", ylab="RMSE", xlab="Percent variance explained by QTL")
par(new=T);plot(t3, c(aut400rmses[2,]), xlim=c(0,10), ylim=c(0,40), type="l", col="blue", ylab="RMSE", xlab="")
par(new=T);plot(t3, c(aut400rmses[3,]), xlim=c(0,10), ylim=c(0,40), type="l", col="green", ylab="RMSE", xlab="")
par(new=T);plot(t3, c(aut400rmses[4,]), xlim=c(0,10), ylim=c(0,40), type="l", col="brown", ylab="RMSE", xlab="")





#####

#arrows(x,yup,xup,y,code=3,angle=90,length=len,col=col,lty=lty,lwd=lwd,...)

save(aut100powersa, aut100powerse, aut100powerss, aut400powersa, aut400powerse, aut400powerss, aut100rmsea, aut100rmsee, aut100rmses, aut400rmsea, aut400rmsee, aut400rmses, file="outpowermse.RData")

save(aut100powersa, aut100powerse, aut100powerss, aut200powersa, aut200powerse, aut200powerss, aut400powersa, aut400powerse, aut400powerss, aut100rmsea, aut100rmsee, aut100rmses, aut200rmsea, aut200rmsee, aut200rmses, aut400rmsea, aut400rmsee, aut400rmses, file="outpowermseN.RData")


load("outpowermse.RData")
par(mfrow=c(2,3))

plot(t1, c(aut100powersa[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="aut100", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powersa[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powersa[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")


plot(1:3, c(aut100powerse[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="equ100", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powerse[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powerse[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")


plot(1:3, c(aut100powerss[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="st100", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powerss[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut100powerss[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")



plot(1:3, c(aut400powersa[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="aut400", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powersa[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powersa[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")


plot(1:3, c(aut400powerse[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="equ400", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powerse[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powerse[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")


plot(1:3, c(aut400powerss[1,]), xlim=c(1,3), ylim=c(0,1), type="l", col="red", main="st400", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powerss[2,]), xlim=c(1,3), ylim=c(0,1), type="l", col="blue", ylab="Power", xlab="")
par(new=T);plot(1:3, c(aut400powerss[3,]), xlim=c(1,3), ylim=c(0,1), type="l", col="green", ylab="Power", xlab="")


#pchisq(0.95, 8)




load("../RDatas/outpowermseN.RData")
x <- list(c(7.992522, 4.417868, 3.085633, 1.504872),
          c(8.391845, 4.371283, 3.111750, 1.569568),
          c(6.775272, 3.618401, 2.013146, 1.175318))

postscript("../Figs/fig5.eps", height=7, width=5, pointsize=12, onefile=FALSE, horizontal=FALSE)

aut100power <- lapply(list(aut100powersa, aut100powerse, aut100powerss), function(a) a*100)
aut200power <- lapply(list(aut200powersa, aut200powerse, aut200powerss), function(a) a*100)
aut400power <- lapply(list(aut400powersa, aut400powerse, aut400powerss), function(a) a*100)

par(mfrow=c(3,3), las=1, xaxs="i", yaxs="i",
    mar=c(5.1, 4.1, 2.1, 1.1), oma=c(0, 3.1, 2.1, 0.0))

for(i in 1:3) {
  plot(x[[i]], aut100power[[i]][4,], xlim=c(0,10), ylim=c(0,100), type="l",
       col="black", main="", ylab="Power",
       xlab="Percent variance explained by QTL")
  for(j in 1:3)
    lines(x[[i]], aut100power[[i]][j,], col=c("red", "blue", "green")[j])

  mtext(side=2, c("Autocorrelated", "Equicorrelated", "Unstructured")[i], line=5.5, las=0, col="darkslateblue")
  if(i==1) mtext(side=3, "n=100", line=2, col="darkslateblue")

  plot(x[[i]], aut200power[[i]][4,], xlim=c(0,10), ylim=c(0,100), type="l",
       col="black", main="", ylab="Power",
       xlab="Percent variance explained by QTL")
  for(j in 1:3)
    lines(x[[i]], aut200power[[i]][j,], col=c("red", "blue", "green")[j])
  if(i==1) mtext(side=3, "n=200", line=2, col="darkslateblue")


  plot(x[[i]], aut400power[[i]][4,], xlim=c(0,10), ylim=c(0,100), type="l",
       col="black", main="", ylab="Power",
       xlab="Percent variance explained by QTL")
  for(j in 1:3)
    lines(x[[i]], aut400power[[i]][j,], col=c("red", "blue", "green")[j])

  if(i==1) mtext(side=3, "n=400", line=2, col="darkslateblue")
}

dev.off()


load("../RDatas/outpowermseN.RData")
x <- list(c(7.992522, 4.417868, 3.085633, 1.504872),
          c(8.391845, 4.371283, 3.111750, 1.569568),
          c(6.775272, 3.618401, 2.013146, 1.175318))

postscript("../Figs/fig6.eps", height=7, width=5, pointsize=12, onefile=FALSE, horizontal=FALSE)

aut100rmse <- list(aut100rmsea, aut100rmsee, aut100rmses)
aut200rmse <- list(aut200rmsea, aut200rmsee, aut200rmses)
aut400rmse <- list(aut400rmsea, aut400rmsee, aut400rmses)

par(mfrow=c(3,3), las=1, xaxs="i", yaxs="i",
    mar=c(5.1, 4.1, 2.1, 1.1), oma=c(0, 3.1, 2.1, 0.0))

for(i in 1:3) {
  plot(x[[i]], aut100rmse[[i]][4,], xlim=c(0,10), ylim=c(0,40),
       type="l", col="black", main="", ylab="RMSE position (cM)", xlab="Percent variance explained by QTL")
  for(j in 1:3)
      lines(x[[i]], aut100rmse[[i]][j,], col=c("red", "blue", "green")[j])

  mtext(side=2, c("Autocorrelated", "Equicorrelated", "Unstructured")[i], line=5.5, las=0, col="darkslateblue")
  if(i==1) mtext(side=3, "n=100", line=2, col="darkslateblue")

  plot(x[[i]], aut200rmse[[i]][4,], xlim=c(0,10), ylim=c(0,40),
       type="l", col="black", main="", ylab="RMSE position (cM)", xlab="Percent variance explained by QTL")
  for(j in 1:3)
      lines(x[[i]], aut200rmse[[i]][j,], col=c("red", "blue", "green")[j])

  if(i==1) mtext(side=3, "n=200", line=2, col="darkslateblue")

  plot(x[[i]], aut400rmse[[i]][4,], xlim=c(0,10), ylim=c(0,40),
       type="l", col="black", main="", ylab="RMSE position (cM)", xlab="Percent variance explained by QTL")
  for(j in 1:3)
      lines(x[[i]], aut400rmse[[i]][j,], col=c("red", "blue", "green")[j])

  if(i==1) mtext(side=3, "n=400", line=2, col="darkslateblue")
}

dev.off()



















