getY <- function(cross, n.max=5, criteria=.9, nn = 0) {

    Y = cross$pheno
    udv <- svd(Y)
    vec <- udv$d^2
    vec <- vec / sum(vec)

    ss = 0

    if( nn == 0) {
        for(j in 1:n.max) {
            ss = ss + vec[j]
            if ( ss > criteria ) {
                nn = j
                break
            }
        }
    }
    pc <- udv$u %*% diag(udv$d)
    pc[,1:nn]
}
