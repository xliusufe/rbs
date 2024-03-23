pvalrp <- function(x,y,rho=0.4,Pk=NULL,criteria=NULL,alpha=0.05,gamma=1.15,family="Fdist",isbic=FALSE){
    n = nrow(x)
    p = ncol(x)
	if(p < rho*n){
        fit = pval(x,y,criteria,alpha,gamma,family,isbic)	
	}
    else{
        k = ceiling(rho*n)
        if(is.null(Pk)){
            Pk = matrix(rnorm((k*p),0,1), nrow = p, ncol = k)
        }
        xp = x%*%Pk
        fit = pval(xp,y,criteria,alpha,gamma,family,isbic)
    }
    fit$rho = rho

    return(fit)
}

