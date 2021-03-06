#' @title Delete or Merge Regressors net
#'
#' @description Fit a path of linear (family = "gaussian") or logistic (family = "binomial") regression models, where models are subsets of continuous predictors and partitions of levels of factors in X. Works even if p>=n (the number of observations is greater than the number of columns in the model matrix).
#'
#' @param X Input data frame; each row is an observation vector; each column can be numerical or integer for a continuous predictor or a factor for a categorical predictor.
#'
#' @param y Response variable; Numerical for family="gaussian" or a factor with two levels for family="binomial". For family="binomial" the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: "gaussian", "binomial".
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details.
#'
#' @param o Parameter of the group lasso screening step, described in Details.
#'
#' @param nlambda Parameter of the group lasso screening step, described in Details.
#'
#' @param lam Value of parameter lambda controling the amount of penalization in rigde regression. Used only for logistic regression in order to allow for parameter estimation in linearly separable setups. Used only for numerical reasons.
#'
#' @param interc Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE). If in X there are any categorical variables, interc=TRUE.
#'
#' @param maxp Maximal number of parameters of the model, smaller values result in quicker computation
#'
#' @details DMRnet algorithm is a generalization of \code{\link{DMR}} to high-dimensional data.
#' It uses a screening step in order to decrease the problem to p<n and DMR subsequently.
#' The screening is done according to the group lasso implemented in the \href{https://CRAN.R-project.org/package=grpreg}{grpreg} package.
#'
#' First, the group lasso for the problem is solved for nlambda values of lambda parameter.
#' Next, for each value of lambda, the scaled nonzero second norms of the groups' coefficients are sorted in decreasing order.
#' Finally, the first i over o fraction of the groups with the largest nonzero values are chosen for further analysis, i = 1,2,...,o-1.
#' E.g., if o=5, first 1/5, first 2/5,..., 4/5 groups with the largest scaled nonzero second norm of coefficients are chosen.
#'
#' The final path of models is chosen by minimizing the likelihood of the models for the number of parameters df equal to 1,...,l<=maxp for some integer l. Note that, in contrast to DMR, the models on the path need not to be nested.
#'
#' @return An object with S3 class "DMR", which  is  a  list  with  the  ingredients:
#'
#' \item{beta}{Matrix p times l of estimated paramters; each column corresponds to a model on the nested path having from l to 1 parameter (denoted as df).}
#' \item{df}{Vector of degrees of freedom; from l to 1.}
#' \item{rss/loglik}{Measure of fit for the nested models: rss (residual sum of squares) for family="gaussian" and loglik (loglikelihood) for family="binomial"}
#' \item{n}{Number of observations.}
#' \item{call}{The call that produced this object.}
#' \item{interc}{If the intercept was fitted: value of parameter interc.}
#'
#' @seealso \code{\link{print.DMR}} for printing, \code{\link{plot.DMR}} for plotting, \code{\link{coef.DMR}} for extracting coefficients and \code{\link{predict.DMR}} for prediction.
#'
#' @examples
#' ## DMRnet for linear regression
#' data(miete)
#' ytr <- miete[1:1500,1]
#' Xtr <- miete[1:1500,-1]
#' Xte <- miete[1501:2053,-1]
#' m1 <- DMRnet(Xtr, ytr)
#' print(m1)
#' plot(m1)
#' g <- gic.DMR(m1, c = 2.5)
#' plot(g)
#' coef(m1, df = g$df.min)
#' ypr <- predict(m1, newx = Xte, df = g$df.min)
#'
#' ## DMRnet for logistic regression
#' data(promoter)
#' ytr <- factor(promoter[1:80,1])
#' Xtr <- promoter[1:80,-1]
#' Xte <- promoter[81:106,-1]
#' m2 <- DMRnet(Xtr, ytr, family = "binomial")
#' print(m2)
#' plot(m2)
#' g <- gic.DMR(m2, c = 2)
#' plot(g)
#' coef(m2, df = g$df.min)
#' ypr <- predict(m2, newx = Xte, df = g$df.min)
#'
#' @export DMRnet

DMRnet <- function(X, y, family = "gaussian", clust.method = "complete", o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){
    X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
    typeofcols <- sapply(1:ncol(X),function(i) class(X[,i]))
    if(sum(unlist(typeofcols) == "ordered") > 0) stop("Error: there is an ordered factor in the data frame, change it to factor")
    sumnonfac <- sum(typeofcols == "factor")
    if (family == "gaussian"){
       if(sumnonfac == 0){
                       return(SOSnet4lm(X, y, o = o, nlambda = nlambda, interc = interc, maxp = maxp))
       } else{
                       return(DMRnet4lm(X, y, clust.method = clust.method, o = o, nlambda = nlambda, maxp = maxp))
       }
    } else{
       if (family == "binomial"){
          if(sumnonfac == 0){
                       return(SOSnet4glm(X, y, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp))
          } else{
                       return(DMRnet4glm(X, y, clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, maxp = maxp))
          }
       }
       else stop("Error: wrong family, should be one of gaussian, binomial")
    }
}
