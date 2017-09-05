#' @title coef.cv.DMR
#'
#' @description Extracts coefficients from a cv.DMR object (for the model with minimal cross-validated error).
#'
#' @param object Fitted cv.DMR object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other coef methods, this function extracts coefficients from a fitted cv.DMR object for the model with minimal cross-validated error.
#'
#' @return Vector of coefficients.
#'
#' @examples
#' ## cv.DMR for linear regression
#' set.seed(13)
#' data(miete)
#' y <- miete$rent
#' X <- miete$area
#' cv = cv.DMR(X,y)
#' coef(cv)
#'
#' @export
coef.cv.DMR <- function(object, ...){
        out <- coef.DMR(object$dmr.fit, df = object$df.min)
        return(out)
}
