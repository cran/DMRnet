DMRnet4glm_help <- function(S, mL, X, y, fl, clust.method, lam){
    if (sum(S) == 0) {
       m <- stats::glm.fit(as.matrix(rep(1, length(y))), y, family = stats::binomial())
       zb = exp(m$coef*rep(1, length(y)))
       pix = zb/(zb + 1)
       loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0])
       b = c(m$coef, sum(fl-1))
       return(list(b = b, loglik = loglik))
    }
    Xn <- X[, S==1, drop = FALSE]
    mfin <- DMR4glm_help(Xn, y, clust.method = clust.method, lam = lam)
    return(mfin)
}
