#' @title gic.DMR
#'
#' @description Computes values of generalized information criterion for the entire sequence of models from a DMR object.
#'
#' @param x Fitted DMR object.
#'
#' @param c Parameter controling amount of penalization for complexity of the model in the generalized information criterion. For linear regression gic for model M is \deqn{GIC_M = RSS_M + df_M*c* log{p}*s^2,} where \eqn{RSS_M} is the residual sum of squares and \eqn{df_M} is the number of parameters in the model M; \eqn{s^2} is an estimator of \eqn{sigma^2} based on the model in the DMR object with the largest number of parameters. For logistic regression gic for model M is \deqn{GIC_M = -2*loglik_M + |M|*c* log{p},} where \eqn{loglik_M} is the logarithm of the likelihood function and \eqn{df_M} is the number of parameters in the model M. Recommended values are c=2.5 for linear regression and c=2 for logistic regression.
#'
#' @return An object of class gic.DMR is returned, which is a list with the ingredients of the gic fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) for the model with minimal gic.}
#'   \item{dmr.fit}{Fitted DMR object.}
#'   \item{gic}{Vector of gic values for the entire sequence of models.}
#' }
#' @seealso \code{\link{plot.gic.DMR}} for plotting, \code{\link{coef.gic.DMR}} for extracting coefficients and \code{\link{predict.gic.DMR}} for prediction.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' (g <- gic.DMR(m, c = 2.5))
#' @export gic.DMR
gic.DMR <- function(x, c = ifelse(x$arguments$family == "gaussian", 2.5, 2)){
        p <- nrow(x$beta)
        n <- x$n
        if (names(x)[3] == "rss"){
           gic <- x$rss + (x$rss[1]/(n - length(x$rss)))*x$df*c*log(p)
        } else{
           gic <- -2*x$loglik + x$df*c*log(p)
        }
        names(gic) <- paste("df", x$df, sep = "")
        kt <- which(gic == min(gic))
        out <- list(df.min = x$df[kt[length(kt)]], dmr.fit = x, gic = gic)
        class(out) <- "gic.DMR"
        return(out)
}
