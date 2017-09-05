w_stats <- function(B, VV, ind1, ind2){
    if (ind1 != ind2){
        n <- ind2 - ind1 + 1
        Tmat <- matrix(0, n, n)
        ind <- t(utils::combn(1:n, 2))
        x <- as.matrix(B[ind[, 1]] - B[ind[, 2]])
        #if (ncol(x) == 1) x <- t(x)
        Tmat[ind] <- (x^2)/(VV[cbind(ind[,1],ind[,1])] + VV[cbind(ind[,2],ind[,2])] - 2*VV[cbind(ind[,1],ind[,2])])
        Tmat <- cbind(0, Tmat)
        t_st <- c(0, B^2/VV[cbind(1:n,1:n)])
        Tmat <- rbind(t_st, Tmat)
    } else {
        Tmat <- matrix(c(0, 0, B^2/VV, 0), 2, 2)
    }
    return(Tmat)
}
