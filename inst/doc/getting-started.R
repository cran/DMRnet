## ---- include = FALSE---------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DMRnet)

## ----miete--------------------------------------------------------------------
data("miete")
X <- miete[,-1]
y <- miete$rent
head(X)

## ----DMRnet-------------------------------------------------------------------
models <- DMRnet(X, y, family="gaussian")
models

## ----plot, fig.height=4, fig.width=6, small.mar=TRUE--------------------------
plot(models, xlim=c(1, 10), lwd=2)

## ----coef---------------------------------------------------------------------
coef(models, df=10)

## ----GIC, fig.height=4, fig.width=6, small.mar=TRUE---------------------------
gic.model <- gic.DMR(models)
plot(gic.model)

## ----gic.df.min---------------------------------------------------------------
gic.model$df.min

## ----gic.coef-----------------------------------------------------------------
coef(gic.model)

## ----cross-validation, fig.height=4, fig.width=6, small.mar=TRUE--------------
cv.model <- cv.DMRnet(X, y)
plot(cv.model)

## ----cv.df.min----------------------------------------------------------------
cv.model$df.min

## ----cv.coef------------------------------------------------------------------
coef(cv.model)==coef(gic.model)

## ----predict------------------------------------------------------------------
predict(gic.model, newx=head(X))
predict(cv.model, newx=head(X))

## ----seq-predict--------------------------------------------------------------
predict(models, newx=head(X))

## ----binomial-----------------------------------------------------------------
binomial_y <- factor(y > mean(y))     #changing Miete response var y into a binomial factor with 2 classes
binomial_models <- DMRnet(X, binomial_y, family="binomial")
gic.binomial_model <- gic.DMR(binomial_models)
gic.binomial_model$df.min

## ----predict-binomial---------------------------------------------------------
predict(gic.binomial_model, newx=head(X), type="class")

