DMR4lm <- function(X, y, clust.method = 'complete'){
    n <- nrow(X)
    if(is.null(colnames(X))) colnames(X) <- paste("x", 1:ncol(X), sep = "")
    if(n != length(y)){
              stop("Error: non-conforming data: nrow(X) not equal to length(y)")
    }
    ssd <- apply(X, 2, function(x) length(unique(x)))
    if (ssd[1] == 1 & (class(X[,1]) == "numeric" | class(X[,1]) == "integer")){
       X <- X[,-1, drop = FALSE]
       ssd <- ssd[-1]
    }
    if(ncol(X) == 0){
              stop("Error: X has zero columns")
    }
    if(sum(ssd == 1) > 0){
               stop("Error: X has columns with sd = 0 apart from the intercept")
    }
    nn <- sapply(1:ncol(X), function(i) class(X[,i]))
    names(nn) <- colnames(X)
    nn[nn == "integer"] <- "numeric"
    if(sum(nn != "numeric" & nn != "factor" ) > 0){
              stop("Error: wrong data type, columns should be one of types: integer, factor, numeric")
    }
    x.full <- stats::model.matrix(y~., data = data.frame(y=y, X[,order(nn), drop = FALSE], check.names = TRUE))
    p <- ncol(x.full)
    if (p >= n){
       stop("Error: p >= n, DMR works only for p < n, use DMRnet instead")
    }
    m <- stats::lm.fit(x.full, y)
    faki <- which(nn == "factor")
    n.factors <- length(faki)
    if (length(faki) > 0){
       n.levels <- sapply(1:n.factors, function(i) length(levels(X[,faki[i]])))
       p.fac <- sum(n.levels - 1)
    } else{
       p.fac <- 0
    }
    cont <- which(nn == "numeric")
    n.cont <- length(cont)
    namCont <- names(nn)[cont]
    ord <- c()
    if(n.cont > 0 ){
              if(n.factors > 0){
                  for (j in 1:n.cont){
                     ord[j] <- sum(n.levels[1:sum(nn[1:cont[j]] == "factor")] - 1) + j + 1
                  }
              } else{
                               ord <- 2:(n.cont + 1)
              }
    }
	  #QR decompostion of the model matrix
    qX <- qr.Q(m$qr)
    rX <- qr.R(m$qr)
    Ro <- solve(rX)
    z <- t(qX)%*%y
    sigma <- as.numeric((t(m$res)%*%m$res)/(n - p))
    #dissimilarity measures - matrices of squared t-statistics for each factor
    if (n.factors > 0){
       Tmats <- lapply(1:n.factors, function(i) {
          i1 <- ifelse(i == 1, 2, sum(n.levels[1:(i - 1)] - 1) + 2)
          i2 <- sum(n.levels[1:i] - 1) + 1
          out <- t_stats(Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = sigma, z = z)
          rownames(out) <- colnames(out) <- m$xlevels[[i]]
          return(out)
       })
       #cutting dendrograms
       models <- lapply(Tmats, function(x) stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL))
       heig <- lapply(1:n.factors, function(x){
            out <- models[[x]]$he
            names(out)<- rep(x, length(out))
            out
       })
       heig <- unlist(heig)
    } else {
        heig <- c()
        models <- list()
    }
    len <- length(heig)
    heig <- c(0,heig)
    names(heig)[1] = "full"
    if ((p.fac + 1) < p){
        if((p.fac + 2) == p){
          heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma*sum(Ro[(p.fac + 2):p,]^2))
        } else {
          heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma*(apply(Ro[(p.fac + 2):p, ], 1, function(y) t(y)%*%y)))
        }
        names(heig.add) <- colnames(x.full)[(p.fac + 2):p]
        heig <- c(heig, heig.add)
    }
    heig <- sort(heig)
    len <- length(heig)
    #fitting models on the path
    Z1 <- Z2 <- c()
    sp <- list()
    form <- c()
    nl <- 0
    Z1 <- X[, faki, drop = FALSE]
    if (n.factors > 0){
            for (i in 1:n.factors){
                sp[[i]] <- 1:n.levels[i]
                sp[[i]][sp[[i]] != 1] <- sp[[i]][sp[[i]] != 1] + nl
                nl <- nl + length(unique(sp[[i]])) - 1
            }
     }
     Z2 <- X[,namCont, drop = FALSE]
     Z <- cbind(Z1,Z2)
     dane <- data.frame(y=y, Z, check.names = T)
     ZZ <- stats::model.matrix(y~., data = dane)
     m <- stats::lm.fit(ZZ, y)
     b <- m$coef
     rss = sum(m$res^2)
     form <- namCont
     if (len > 2){
       for (i in 2:(len - 1)){
         kt <- names(heig)[i]
         if(length(intersect(kt, namCont)) > 0){
                          form <- form[-which(form == kt)]
                          Z2 <- Z2[, form, drop = FALSE]

         } else {
           kt <- as.numeric(kt)
           dod <- min(sp[[kt]][sp[[kt]] != 1])
           sp[[kt]] <- stats::cutree(models[[kt]], h = heig[i])
           if(length(sp[[kt]][sp[[kt]] != 1]) > 0){
                                       sp[[kt]][sp[[kt]] != 1] <- sp[[kt]][sp[[kt]] != 1] + dod - min(sp[[kt]][sp[[kt]] != 1])
           }
           Z1[,kt] <- X[, faki[kt]]
           levels(Z1[,kt]) <- sp[[kt]]
           Z1[,kt] <- factor(Z1[,kt])
           if (kt < length(sp)) for( x in (kt+1):length(sp)){ if (length(sp[[x]][sp[[x]]!=1]) > 0 ) sp[[x]][sp[[x]]!= 1] = sp[[x]][sp[[x]]!=1] - 1}
           nl <- nl - 1
         }
         Z <- cbind(Z1[,which(apply(Z1, 2, function(x) length(unique(x))) != 1)], Z2)
         dane <- data.frame(y = y, Z, check.names = T)
         ZZ <- stats::model.matrix(y~., data = dane)
         m <- stats::lm.fit(ZZ, y)
         be <- c(0, m$coef[-1])
         bb <- c()
         if(n.factors > 0){
                       bb <- unlist(sapply(1:length(sp), function(j) sp[[j]][-1]))
         }
         bb2 <- rep(1, n.cont)
         names(bb2) <- namCont
         if(length(form) > 0){
                    bb2[form] <- (nl + 2):(nl + 1 + length(form))
         }
         bb <- c(bb, bb2)
         b=cbind(b, c(m$coef[1], be[bb]))
         rss = c(rss, sum(m$res^2))
     }
   }
   m <- stats::lm.fit(as.matrix(rep(1, length(y))), y)
   b = cbind(b, c(m$coef, rep(0, length(heig) - 1)))
   rss = c(rss, sum(m$res^2))
   if(length(ord) > 0){
                  ind <- c(1:p)
                  ind[-ord] = 1:(p - length(ord))
                  ind[ord] = (p - length(ord) + 1):p
                  b = b[ind,]
   }
   fit <- list(beta = b, df = p:1, rss = rss, n = n, arguments = list(family = "binomial", clust.method = clust.method), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
