#' @title cross-validation for DMRnet
#'
#' @description Does k-fold cross-validation for DMR and returns a value for df.
#'
#' @param X Input data frame, of dimension n x p; each row is an observation vector. Columns can be numerical or integer for continuous predictors or factors for categorical predictors.
#'
#' @param y Response variable. Numerical for family="gaussian" or a factor with two levels for family="binomial". For family="binomial" the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: "gaussian", "binomial".
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details.
#'
#' @param o Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param nlambda Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param lam Value of parameter lambda controling the amount of penalization in rigde regression. Used only for logistic regression in order to allow for parameter estimation in linearly separable setups. Used only for numerical reasons.
#'
#' @param interc Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE). If in X there are any categorical variables, interc=TRUE.
#'
#' @param nfolds Number of folds in cross-validation.
#'
#' @param maxp Maximal number of parameters of the model, smaller values result in quicker computation.
#'
#' @details cv.DMRnet algorithm does k-fold cross-validation for DMRnet. The df for the minimal estimated prediction error is returned.
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
#' ## cv.DMRnet for linear regression
#' set.seed(13)
#' data(miete)
#' ytr <- miete$rent[1:1500]
#' Xtr <- miete$area[1:1500]
#' Xte <- miete$area[1501:2053]
#' cv <- cv.DMRnet(Xtr, ytr)
#' print(cv)
#' plot(cv)
#' coef(cv)
#' ypr <- predict(cv, newx = Xte)
#'
#' @export cv.DMRnet

cv.DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){
       X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
       if (family == "gaussian"){
          n <- length(y)
          foldid <- cvfolds(n, nfolds)
          error <- list()
          for (fold in 1:nfolds){
              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMRnet(Xtr, ytr, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
              pred <- predict.DMR(dmr, newx = as.data.frame(Xte))
              error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
          }
          foldmin <- min(sapply(error, length))
          error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
          error <- rowSums(error)/n
          dmr.fit <- DMRnet(X, y, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
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
          error <- list()
          for (fold in 1:nfolds){
              Xte <- X[foldid == fold, , drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, , drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMRnet(Xtr, ytr, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
              pred <- predict.DMR(dmr, newx = as.data.frame(Xte), type = "class")
              error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
          }
          foldmin <- min(sapply(error, length))
          error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
          error <- rowSums(error)/(n1 + n2)
          dmr.fit <- DMRnet(X, y, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
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
