cvfolds <- function (n, K = 5){
    if (!(n > 0))
        stop("'n' must be positive")
    if (!((K > 1) && K <= n))
        stop("'K' outside allowable range")
    subsets <- sample(n)
    which <- rep(seq_len(K), length.out = n)
    block <- which[order(subsets)]
    return(block)
}