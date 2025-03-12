sirs <- function(Y, X, standardize_X = TRUE, N = NULL, d = 10, ntop = 10){
    if (standardize_X){
        X = scale(X);
    }

    n = nrow(X)
    p = ncol(X)
    if(is.null(n)) n = length(X)
    if(is.null(p)) p = 1

    if (is.null(N)) {
        N = floor(n / log(n))
    }
    if (N > p) N = p;

    dims = c(n,p)
    fit <- .Call("_SIRS", 
                as.numeric(Y),
                as.numeric(t(X)), 
                as.integer(dims), 
                as.integer(N),
                as.integer(d),
                as.integer(ntop)
            )

    fit$indn <- fit$indn + 1
    return(fit)
}
