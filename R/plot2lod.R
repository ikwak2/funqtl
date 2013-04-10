plot2lod <-
function(cross, z2, ylab="QTL position", xlab="Time", mval=0, ...) {
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
    image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab, ...)

    last <- end[nch]
    start <- c(start, last + 2)
    for(i in 1:nch)
        rect(0,end[i], nrow(z2)+1,start[i+1], col="white")

    rect(0,end[i], nrow(z2)+1,start[i+1]+1, col="white")


    u <- par("usr") # plot ranges [left,right, bottom,top]

    width <- 0.01*diff(u[1:2])
    text(u[1]-5*width, midpt, 1:nch, xpd=TRUE)

    for(i in seq(along=start)) {
        if(i %% 2)
            rect(u[1]-width, start[i], u[1]-width*2, end[i], col="gray30", xpd=TRUE)
        else
            rect(u[1]-width*2, start[i], u[1]-width*3, end[i], col="gray30", xpd=TRUE)
    }
}
