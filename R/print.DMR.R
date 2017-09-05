#' @title print.DMR
#'
#' @description Print a DMR object.
#'
#' @param x Fitted DMR object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Print a summary of the DMR path at each step along the path.
#'
#' @return The summary is silently returned.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' print(m)
#' @export
print.DMR <- function (x,...){
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (names(x)[3] == "rss"){
       print(cbind(Df = x$df, RSS = x$rss))
    } else {
       print(cbind(Df = x$df, loglik = x$loglik))
    }
}
