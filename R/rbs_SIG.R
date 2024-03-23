rbs_sig <- function(x,y,V=NULL,gamma=1.5, lambda=NULL,criteria=2,nflip=1,tau=1){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
    
    isV = ifelse(is.null(V),0,1)
    if(isV==0) V = diag(q)
    
    if(criteria==0){
        if(is.null(lambda)){
            alpha0 = 0.05/q
            lambda = qchisq(alpha0,p,lower.tail = F)
        } 
        dims    = c(n,p,q,isV,nflip)
        param   = c(gamma,lambda)
        fit     <- .Call("RBSS_FLIP", 
                        as.numeric(x), 
                        as.numeric(y), 
                        as.numeric(V), 
                        as.integer(dims), 
                        as.numeric(param)
                    )
    }
    else{
        if(is.null(lambda)){
            alpha0  = 0.05/q
            threh   = qchisq(alpha0,p,lower.tail = F)
            lambda  =  seq(threh/5, threh*5, length = 50)  
        } 
        nlam    = length(lambda)
        dims    = c(n,p,q,nlam,isV,nflip,criteria)
        param   = c(gamma,tau)
        fit     <- .Call("RBSS_FLIP_BIC", 
                        as.numeric(x), 
                        as.numeric(y), 
                        as.numeric(V), 
                        as.numeric(lambda), 
                        as.integer(dims), 
                        as.numeric(param)
                    )
        selected        <- fit$selected + 1
        fit$selected    <- selected
        deltahat        <- matrix(fit$delta,q,nlam)
        fit$deltapath   <- deltahat
        fit$delta       <- deltahat[,selected]       
    }
    
    return(fit)
}