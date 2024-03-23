rbs <- function(x,y,gamma=1.5, lambda=NULL,criteria=2,tau=1){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
     
    if(criteria==0){
        if(is.null(lambda)){
            alpha0 = 0.05/q
            lambda = qchisq(alpha0,p,lower.tail = T)
        } 
        dims = c(n,p,q)
        param = c(gamma,lambda)
        fit <- .Call("RBSS", 
                    as.numeric(x), 
                    as.numeric(y), 
                    as.integer(dims), 
                    as.numeric(param)
                )
    }
    else{
        if(is.null(lambda)){
            alpha0 = 0.05/q
            threh = qchisq(alpha0,p,lower.tail = T)
            lambda =  seq(threh/5, threh*5, length = 50)  
        } 
        nlam = length(lambda)
        dims = c(n,p,q,nlam,criteria)
        param = c(gamma,tau)
        fit <- .Call("RBSS_BIC", 
                    as.numeric(x), 
                    as.numeric(y), 
                    as.numeric(lambda), 
                    as.integer(dims), 
                    as.numeric(param)
                )
        selected <- fit$selected + 1
        deltahat <- matrix(fit$delta,q,nlam)[,selected]
        fit$deltapath <- fit$delta
        fit$delta <- deltahat
        fit$selected <- selected
    }

    return(fit)
}
