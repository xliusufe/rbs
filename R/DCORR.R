dcorr <- function(y, x){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
    
    dims = c(n,q,p)
    fit <- .Call("DCORR", y, x, as.integer(dims))
    
    fit
}

sisdc <- function(y, x, d=1, ntop=10){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
    
    if(d==1){
        dims = c(n,q,p)
        fit <- .Call("DCORR_SCREEN", y, x, as.integer(dims), as.integer(ntop))
    }
    else if(d==2){
        dims = c(n,p,q)
        fit <- .Call("DCORR_SCREEN", x, y, as.integer(dims), as.integer(ntop))
    }
    else print("d must be 1 or 2!")
    
    fit$indn <- fit$indn + 1
    fit
}