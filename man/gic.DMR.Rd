% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gic.DMR.R
\name{gic.DMR}
\alias{gic.DMR}
\title{gic.DMR}
\usage{
gic.DMR(x, c = ifelse(x$arguments$family == "gaussian", 2.5, 2))
}
\arguments{
\item{x}{Fitted DMR object.}

\item{c}{Parameter controling amount of penalization for complexity of the model in the generalized information criterion. For linear regression gic for model M is \deqn{GIC_M = RSS_M + df_M*c* log{p}*s^2,} where \eqn{RSS_M} is the residual sum of squares and \eqn{df_M} is the number of parameters in the model M; \eqn{s^2} is an estimator of \eqn{sigma^2} based on the model in the DMR object with the largest number of parameters. For logistic regression gic for model M is \deqn{GIC_M = -2*loglik_M + |M|*c* log{p},} where \eqn{loglik_M} is the logarithm of the likelihood function and \eqn{df_M} is the number of parameters in the model M. Recommended values are c=2.5 for linear regression and c=2 for logistic regression.}
}
\value{
An object of class gic.DMR is returned, which is a list with the ingredients of the gic fit.
\describe{
  \item{df.min}{df (number of parameters) for the model with minimal gic.}
  \item{dmr.fit}{Fitted DMR object.}
  \item{gic}{Vector of gic values for the entire sequence of models.}
}
}
\description{
Computes values of generalized information criterion for the entire sequence of models from a DMR object.
}
\examples{
data(miete)
y <- miete[,1]
X <- miete[,-1]
m <- DMR(X, y)
(g <- gic.DMR(m, c = 2.5))
}
\seealso{
\code{\link{plot.gic.DMR}} for plotting, \code{\link{coef.gic.DMR}} for extracting coefficients and \code{\link{predict.gic.DMR}} for prediction.
}

