rbs_sig <- function(x,y,V=NULL,gamma=1.5, lambda=NULL,criteria=2,nflip=1,tau=1,is_Bdiag=TRUE){
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
            lambda = qchisq(alpha0,p,lower.tail = T)
        } 
        dims = c(n,p,q,isV,nflip)
        param = c(gamma,lambda)
        if(is_Bdiag){
            fit <- .Call("RBSS_FLIP", x, y, V, as.integer(dims), param)    
        }else{
            fit <- .Call("RBSS_FLIP_ND", x, y, V, as.integer(dims), param)
        }
        
        
    }
    else{
        if(is.null(lambda)){
            alpha0 = 0.05/q
            threh = qchisq(alpha0,p,lower.tail = T)
            lambda =  seq(threh/5, threh*5, length = 50)  
        } 
        nlam = length(lambda)
        dims = c(n,p,q,nlam,isV,nflip,criteria)
        param = c(gamma,tau)
        if(is_Bdiag){
            fit <- .Call("RBSS_FLIP_BIC", x, y, V, lambda, as.integer(dims), param)
        }else{
            fit <- .Call("RBSS_FLIP_ND_BIC", x, y, V, lambda, as.integer(dims), param)
        }
        selected <- fit$selected + 1
        fit$selected <- selected
        deltahat <- matrix(fit$delta,q,nlam)[,selected]
        fit$deltapath <- fit$delta
        fit$delta <- deltahat        
    }
    
    fit
}