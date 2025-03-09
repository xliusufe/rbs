pcorr <- function(Y, X){
    n = nrow(X)
    p = ncol(X)
    if(is.null(n)) n = length(X)
    if(is.null(p)) p = 1

    dims = c(n,p)
    fit <- .Call("PEARSON_CORR", 
                as.numeric(Y),
                as.numeric(t(X)),
                as.integer(dims)
            )
    
    return(fit)
}

sis <- function(Y, X, ntop=10){
    n = nrow(X)
    p = ncol(X)
    if(is.null(n)) n = length(X)
    if(is.null(p)) p = 1
    
    dims = c(n,p)
    fit <- .Call("PEARSON_SCREEN", 
                as.numeric(Y),
                as.numeric(t(X)),
                as.integer(dims), 
                as.integer(ntop)
            )

    fit$indn <- fit$indn + 1
    return(fit)
}