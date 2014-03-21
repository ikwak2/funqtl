#' Plot LOD matrix.
#'
#' Plot a matrix of LOD scores, for which the (i,j)th element is for phenotype
#' i at marker j.
#'
#'
#' @param cross An object of class 'cross'. See 'read.cross' for details.
#' @param out1 A scanone output of LOD scores.
#' @param ylab A lable of y axis. Default is 'QTL position'.
#' @param xlab A lable of x axis. Default is 'Time'.
#' @param mval The maximum LOD value of legend. The color of legend goes 0 to
#' 'mval'. If this value is less than the max lod score, it automatically
#' changed to max value.
#' @param col Vector of colors
#' @param \dots More graphical components of 'image.plot'.
#' @return A graph of LOD matrix.
#' @author Il-Youp Kwak, <email: ikwak2@@stat.wisc.edu>
#' @seealso 'plotlodmatlist', 'plotlodmatlist2'
#' @keywords hplot
#' @examples
#'
#'
#' # call dataset
#' data(simspal)
#' simspal <- calc.genoprob(simspal, step=1)
#' out1 <- scanone(simspal, pheno.col=1:241 , method="hk")
#'
#' \dontshow{par(mar=c(3.1, 3.1, 1.6, 0.6))}
#' plot2lod(simspal, out1, main="The LOD image of the data")
#'
#'

plot2lod <-
function (cross, out1, ylab = "QTL position", xlab = "Time", mval = 0,
    col = rev(heat.colors(100)), ...)
{
    z2 <- t(as.matrix(out1[,c(-(1:2))]))
    nch <- nchr(cross)
    nchs <- table(out1[,1])

    end <- NULL
    start <- NULL
    num = 0
    z3 <- NULL
    for (i in 1:nch) {
        z3 <- cbind(z3, cbind(z2[, ((num + 1):(num + nchs[i]))],
            rep(mval, nrow(z2)), rep(mval, nrow(z2))))
        num = num + nchs[i]
        end <- c(end, num + (i - 1) * 2)
        start <- c(start, num - nchs[i] + 1 + (i - 1) * 2)
    }

    midpt <- (start + end)/2
    x <- 1:(dim(z3)[1])
    y <- 1:(dim(z3)[2])

    image.plot(x, y, z3, yaxt = "n", xlab = xlab, ylab = ylab,
        col = col, ...)
    last <- end[nch]
    start <- c(start, last + 2)
    for (i in 1:nch) rect(0, end[i], nrow(z2) + 1, start[i +
        1], col = "white")
    rect(0, end[i], nrow(z2) + 1, start[i + 1] + 1, col = "white")
    u <- par("usr")
    width <- 0.01 * diff(u[1:2])
    text(u[1] - 5 * width, midpt, 1:nch, xpd = TRUE)
    for (i in seq(along = start)) {
        if (i%%2)
            rect(u[1] - width, start[i], u[1] - width * 2, end[i],
                col = "gray30", xpd = TRUE)
        else rect(u[1] - width * 2, start[i], u[1] - width *
            3, end[i], col = "gray30", xpd = TRUE)
    }
}
