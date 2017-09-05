#' @title cross-validation for DMR
#'
#' @description Does k-fold cross-validation for DMR and returns a value for df.
#'
#' @param X Input data frame, of dimension n x p; DMR works only if p<n, for p>=n see DMRnet; each row is an observation vector. Columns can be numerical or integer for continuous predictors or factors for categorical predictors.
#'
#' @param y Response variable. Numerical for family="gaussian" or a factor with two levels for family="binomial". For family="binomial" the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: "gaussian", "binomial".
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details.
#'
#' @param lam Value of parameter lambda controling the amount of penalization in rigde regression. Used only for logistic regression in order to allow for parameter estimation in linearly separable setups.
#'
#' @param nfolds Number of folds in cross-validation.
#'
#' @details cv.DMR algorithm does k-fold cross-validation for DMR. The df for the minimal estimated prediction error is returned.
#'
#' @return An object with S3 class "cv.DMR" is  returned,  which  is  a  list  with  the  ingredients  of  the  cross-validation fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) for the model with minimal cross-validated error.}
#'   \item{dmr.fit}{Fitted DMR object for the full data.}
#'   \item{cvm}{The mean cross-validated error for the entire sequence of models.}
#'   \item{foldid}{The fold assignments used.}
#' }
#'
#' @seealso  \code{\link{plot.cv.DMR}} for plotting, \code{\link{coef.cv.DMR}} for extracting coefficients and \code{\link{predict.cv.DMR}} for prediction.
#'
#' @examples
#' ## cv.DMR for linear regression
#' set.seed(13)
#' data(miete)
#' ytr <- miete$rent[1:1500]
#' Xtr <- miete$area[1:1500]
#' Xte <- miete$area[1501:2053]
#' cv <- cv.DMR(Xtr, ytr)
#' print(cv)
#' plot(cv)
#' coef(cv)
#' ypr <- predict(cv, newx = Xte)
#'
#' @export cv.DMR

cv.DMR <- function(X, y, family = "gaussian", clust.method = 'complete', lam = 10^(-7), nfolds = 10){
       X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
       if (family == "gaussian"){
          n <- length(y)
          foldid <- cvfolds(n, nfolds)
          error <- c()
          for (fold in 1:nfolds){
              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMR(Xtr, ytr, family = "gaussian", clust.method = clust.method)
              pred <- predict.DMR(dmr, newx = as.data.frame(Xte))
              error <- cbind(error, apply(pred, 2, function(z) sum((z - yte)^2)))
          }
          error <- rowSums(error)/n
          dmr.fit <- DMR(X, y, family = "gaussian", clust.method = clust.method)
          kt <- which(error == min(error))
          df.min <- dmr$df[kt[length(kt)]]
       } else{
         if (family == "binomial"){
          if (class(y) != "factor"){
             stop("Error: y should be a factor")
          }
          lev <- levels(factor(y))
          if (length(lev) != 2){
             stop("Error: factor y should have 2 levels")
          }
          n1 <- table(y)[1]
          n2 <- table(y)[2]
          foldid1 <- cvfolds(n1, nfolds)
          foldid2 <- cvfolds(n2, nfolds)
          foldid <- c()
          foldid[which(y == levels(factor(y))[1])] = foldid1
          foldid[which(y == levels(factor(y))[2])] = foldid2
          error <- c()
          for (fold in 1:nfolds){
              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMR(Xtr, ytr, family = "binomial", clust.method = clust.method, lam = lam)
              pred <- predict.DMR(dmr, newx = as.data.frame(Xte), type = "class")
              error <- cbind(error, apply(pred, 2, function(z) sum(z != yte)))  # moze sie wysypac!!! jezeli roznej wielkosci dfy
          }
          error <- rowSums(error)/(n1 + n2)
          dmr.fit <- DMR(X, y, family = "binomial", clust.method = clust.method, lam = lam)
          kt <- which(error == min(error))
          df.min <- dmr$df[kt[length(kt)]]
         }
         else{
              stop("Error: wrong family, should be one of gaussian, binomial")
         }
       }
       out <- list(df.min = df.min, dmr.fit = dmr.fit, cvm = error, foldid = foldid)
       class(out) <- "cv.DMR"
       return(out)
}
