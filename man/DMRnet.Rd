% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMRnet.R
\name{DMRnet}
\alias{DMRnet}
\title{Delete or Merge Regressors net}
\usage{
DMRnet(
  X,
  y,
  family = "gaussian",
  o = 5,
  nlambda = 100,
  lam = 10^(-7),
  interc = TRUE,
  maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4)),
  lambda = NULL,
  algorithm = "DMRnet",
  clust.method = ifelse(algorithm == "glamer", "single", "complete")
)
}
\arguments{
\item{X}{Input data frame; each row is an observation vector; each column can be numerical or integer for a continuous predictor or a factor for a categorical predictor.}

\item{y}{Response variable; Numerical for \code{family="gaussian"} or a factor with two levels for \code{family="binomial"}. For \code{family="binomial"} the last level in alphabetical order is the target class.}

\item{family}{Response type; one of: \code{"gaussian"}, \code{"binomial"}.}

\item{o}{Parameter of the group lasso screening step, described in Details, the default value is 5.}

\item{nlambda}{Parameter of the group lasso screening step, described in Details, the default value is 100.}

\item{lam}{The amount of penalization in ridge regression (used for logistic regression in order to allow for parameter estimation in linearly separable setups) or the amount of matrix regularization in case of linear regression. Used only for numerical reasons. The default is 1e-7.}

\item{interc}{Should intercept(s) be fitted (the default, \code{interc=TRUE}) or set to zero (\code{interc=FALSE}). If in \code{X} there are any categorical variables, \code{interc=TRUE} must be set.}

\item{maxp}{Maximal number of parameters of the model, smaller values result in quicker computation}

\item{lambda}{Explicitly provided net of lambda values for the group lasso screening step, described in Details. If provided, it overrides the value of \code{nlambda} parameter.}

\item{algorithm}{The algorithm to be used; for partition selection (merging levels) use one of: \code{"DMRnet"} (the default), \code{"glamer"} or \code{"PDMR"}. Alternatively, use \code{"var_sel"} for variable (group) selection with no partition selection.}

\item{clust.method}{Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details. \code{clust.method="complete"} is the default for all algorithms except \code{algorithm="glamer"}, for which \code{clust.method="single"} is the default.}
}
\value{
An object with S3 class \code{"DMR"}, which  is  a  list  with  the  ingredients:

\item{beta}{Matrix p times l of estimated parameters; each column corresponds to a model on the nested path having from l to 1 parameter (denoted as df).}
\item{df}{Vector of degrees of freedom; from l to 1.}
\item{rss/loglik}{Measure of fit for the nested models: rss (residual sum of squares) is returned for \code{family="gaussian"} and loglik (loglikelihood) is returned for \code{family="binomial"}.}
\item{n}{Number of observations.}
\item{levels.listed}{Minimal set of levels of respective factors present in data.}
\item{lambda}{The net of lambda values used in the screening step.}
\item{arguments}{List of the chosen arguments from the function call.}
\item{interc}{If the intercept was fitted: value of parameter \code{interc} is returned.}
}
\description{
Fits a path of linear (\code{family="gaussian"}) or logistic (\code{family="binomial"}) regression models, where models are subsets of continuous predictors and partitions of levels of factors in \code{X}. Works even if p>=n (the number of observations is greater than the number of columns in the model matrix).
}
\details{
\code{DMRnet} algorithm is a generalization of \code{\link{DMR}} to high-dimensional data.
It uses a screening step in order to decrease the problem to p<n and then uses \code{DMR} subsequently.
The screening is done with the group lasso algorithm implemented in the \href{https://CRAN.R-project.org/package=grpreg}{grpreg} package.

First, the group lasso for the problem is solved for \code{nlambda} values of lambda parameter, or for the net of lambda values (if \code{lambda} is explicitly provided).
Next, for each value of lambda, the scaled nonzero second norms of the groups' coefficients are sorted in decreasing order.
Finally, the first i over \code{o} fraction of the groups with the largest nonzero values are chosen for further analysis, i = 1,2,...,\code{o}-1.
E.g., if \code{o}=5, first 1/5, first 2/5,..., 4/5 groups with the largest scaled nonzero second norm of coefficients are chosen.

The final path of models is chosen by minimizing the likelihood of the models for the number of parameters df equal to 1,...,l<=\code{maxp} for some integer l. Note that, in contrast to \code{DMR}, the models on the path need not to be nested.
}
\examples{
## DMRnet for linear regression
data(miete)
ytr <- miete[1:200,1]
Xtr <- miete[1:200,-1]
Xte <- miete[201:250,-1]
m1 <- DMRnet(Xtr, ytr)
print(m1)
plot(m1)
g <- gic.DMR(m1, c = 2.5)
plot(g)
coef(m1, df = g$df.min)
ypr <- predict(m1, newx = Xte, df = g$df.min)

## DMRnet for logistic regression
data(promoter)
ytr <- factor(promoter[1:70,1])
Xtr <- promoter[1:70,-1]
Xte <- promoter[71:106,-1]
m2 <- DMRnet(Xtr, ytr, family = "binomial")
print(m2)
plot(m2)
g <- gic.DMR(m2, c = 2)
plot(g)
coef(m2, df = g$df.min)
ypr <- predict(m2, newx = Xte, df = g$df.min)

## PDMR for linear regression
data(miete)
ytr <- miete[1:200,1]
Xtr <- miete[1:200,-1]
Xte <- miete[201:250,-1]
m1 <- DMRnet(Xtr, ytr, algorithm="PDMR")
print(m1)
plot(m1)
g <- gic.DMR(m1, c = 2.5)
plot(g)
coef(m1, df = g$df.min)
ypr <- predict(m1, newx = Xte, df = g$df.min)

}
\seealso{
\code{\link{print.DMR}} for printing, \code{\link{plot.DMR}} for plotting, \code{\link{coef.DMR}} for extracting coefficients and \code{\link{predict.DMR}} for prediction.
}
