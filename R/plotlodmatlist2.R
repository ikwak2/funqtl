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