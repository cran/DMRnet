#' @title coef.DMR
#'
#' @description Extracts coefficients from a DMR object.
#'
#' @param object Fitted DMR object.
#'
#' @param df Number of parameters in the model for which coefficients are required. Default is the entire path of models.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other coef methods, this function extracts coefficients from a fitted DMR object.
#'
#' @return Vector or matrix of coefficients.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' coef(m, df = 12)
#' @export
coef.DMR <- function(object, df = NULL, ...){
         if(is.null(df)){
                         out <- object$beta
                         colnames(out) <- paste("df", object$df, sep = "")
                         return(out)
         }
         if(class(df) != "numeric" & class(df) != "integer"){
                      stop("Error: wrong input type, df should have type numeric")
         }
         out <- object$beta[,ncol(object$beta) - df + 1]
         return(out)
}
