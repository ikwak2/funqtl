getprofile <-
function (tpy=c("comb","sep"), ...) {

    tpy <- match.arg(tpy)

    if( tpy == "comb" ) {
        out <- profileLodMatfn2(...)
    } else if (tpy == "sep") {
        out <- profileLodMatfn(...)
    }
    out
}
