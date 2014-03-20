plotprofile <- function(lodmatlist, ylab="QTL position", xlab="Time", mval=0, ...) {

    if(class(lodmatlist) == "lodprofileM") {
        plotlodmatlist(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    } else if (class(lodmatlist) == "lodprofileM2") {
        plotlodmatlist2(lodmatlist, ylab=ylab, xlab=xlab, mval=mval, ...)
    }
}
