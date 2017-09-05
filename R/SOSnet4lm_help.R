SOSnet4lm_help <- function(S, mL, X, y, interc = interc){
    screenPred <- which(S == 1)
    s <- sum(S)
    if (interc == FALSE){
       QR <- qr(X[, screenPred, drop = FALSE])
    } else{
       X <- cbind(1, X)
       QR <- qr(cbind(1, X[, screenPred + 1, drop = FALSE]))
       s <- s + 1
    }
    W <- qr.R(QR)
    Q <- qr.Q(QR)
    z <- t(Q)%*%y
    bety <- backsolve(W, z)
    W_ <- backsolve(W, diag(s))
    v_bety <- apply(W_^2, 1, sum)
    T2 <- bety^2/v_bety
    if (interc == TRUE) T2 <- T2[-1]
    ind <- sort(T2, decreasing = TRUE, index.return = TRUE)$ix
    if (interc == TRUE){
       QR <- qr(X[, c(1, screenPred[ind] + 1), drop = FALSE])
    } else{
       QR <- qr(X[, screenPred[ind], drop = FALSE])
    }
    q2 <- (t(qr.Q(QR))%*%y)^2
    rss <- rep(NA, s)
    rss1 <- sum(y^2) - q2[1]
    rss[s] <- rss1
    if (s > 1){
       for (k in (s-1):1){
           rss[k] <- rss[k+1] - q2[s-k+1]
       }
    }
    return(list(rss = rss, ind = screenPred[ind]))
}
