#' @title plot.gic.DMR
#'
#' @description Plot gic values from a gic.DMR object.
#'
#' @param x Fitted gic.DMR object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a plot of generalized information criterion for the entire sequence of models from the fitted gic.DMR object.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' g <- gic.DMR(m, c = 2.5)
#' plot(g)
#' @export
plot.gic.DMR <- function(x, ...){
  graphics::plot(x$dmr.fit$df, x$gic, pch = 16, col = 1, xlab = "df", ylab = "GIC")
  graphics::points(x$df.min, min(x$gic), pch = 16, col = 2)
}
