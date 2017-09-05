#' @title plot.DMR
#'
#' @description Plot coefficients from a DMR object.
#'
#' @param x Fitted DMR object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a coefficient profile plot of the coefficient paths for a fitted DMR object.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' plot(m)
#' @export
plot.DMR <- function (x, ...){
         if (x$interc == TRUE){
           graphics::matplot(x$df, t(x$beta[-1,]), type = "l", lty = "solid", ylab = "Coefficients", xlab = "df")
         } else{
           graphics::matplot(x$df, t(x$beta), type = "l", lty = "solid", ylab = "Coefficients", xlab = "df")
         }
}
