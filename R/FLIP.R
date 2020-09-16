flip <- function(A,b=NULL,x0=NULL,nflip=1){
    q = ncol(A)
    if(is.null(x0)) x0 = rep(0,q)
    if(is.null(b)) b = rep(0.0,q)
    param = c(q,nflip)
    fit <- .Call("FLIP", as.numeric(A), b, as.integer(x0), as.integer(param))
    
    fit
}
