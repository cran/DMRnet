DMRnet4glm <- function(X, y, clust.method, o, nlambda, lam, maxp, lambda) {

    y <- prelasso_binomial(y)

    out <- prelasso_common(X, y)
    n <-             out$n
    levels.listed <- out$levels.listed
    fl <-            out$fl
    x.full <-        out$x.full_normalized
    p <-             out$p
    p.x <-           out$p.x
    ord <-           out$ord
    groups <-        out$groups
    X <-             out$X

    if (is.null(lambda)) {
      user.lambda<-substitute()    #make user.lambda - paradoxically - not present in a call to grpreg
    } else {
      nlambda <- length(lambda)   #override this parameter
      user.lambda <- lambda
    }

    mL <-  grpreg::grpreg(x.full[,-1, drop=FALSE], y, group=groups, penalty = "grLasso", family ="binomial", nlambda = nlambda, lambda = user.lambda)

    bb <-  postlasso_common(mL$lambda, n, mL$beta)
    fac <- postlasso_fac(bb, groups)   #fac must be computed on bb without intercept, it happens internally in postlasso_fac()
    SS <-  postlasso_O_step_preparation(p, p.x, n, o, fac, interc=TRUE)

    mm <-  lapply(1:ncol(SS), function(i) DMRnet4glm_help(SS[,i], X, y, fl, clust.method = clust.method, lam = lam))

    return(wrap_up_binomial(mm, p, maxp, SS, fl, x.full, ord, n, levels.listed, mL, list(family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, maxp = maxp, lambda = lambda)))
}
