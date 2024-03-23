rbsrp <- function(x,y,rho=0.4,Pk=NULL,gamma=1.5, lambda=NULL,criteria=2,tau=1){
    n = nrow(x)
    p = ncol(x)
	if(p < rho*n){
        fit = rbs(x,y,gamma,lambda,criteria,tau)	
	}
    else{
        k = ceiling(rho*n)
        if(is.null(Pk)){
            Pk = matrix(rnorm((k*p),0,1), nrow = p, ncol = k)
        }
        xp = x%*%Pk
        fit = rbs(xp,y,gamma,lambda,criteria,tau)
    }
    fit$rho = rho

    return(fit)
}
