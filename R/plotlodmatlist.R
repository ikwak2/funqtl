#' @importFrom graphics par rect text

plotlodmatlist <- function(lodmatlist, times, ylab="QTL position", xlab="Time", mval=0, ...) {


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
    if(!missing(times) && !is.null(times)) {
      if(length(times) != length(x))
        stop("times should have length ", length(x))
      x <- times
    }

    par(mar=c(5.1,6.1,2.1,2.1))
    image.plot(x, y, z3, yaxt="n", xlab=xlab, ylab=ylab, ...)

    last <- end[nlst]
    start <- c(start, last + 2)
    for(i in 1:nlst)
      rect(x[1]-diff(x[1:2]),end[i], x[length(x)]+diff(x[1:2]),start[i+1]-0.5, col="white", border="white")
    rect(x[1]-diff(x[1:2]),max(end), x[length(x)]+diff(x[1:2]),max(start)+1, col="white", border="white")

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
