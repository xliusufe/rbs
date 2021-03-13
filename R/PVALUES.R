pval <- function(x,y,criteria=NULL,alpha=0.05,gamma=1.15,family="Fdist",isbic=FALSE){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
    
    ngamma = length(gamma)
    dims = c(n,p,q)
    fit = .Call("PVALUES",x,y,as.integer(dims))
    if(family=="Chi2"){
        if(is.null(criteria) || criteria=="RBS"){
            if(isbic && ngamma>1){
                dims    = c(n,p,q,ngamma)
                fit     = .Call("PVALUES_BIC",x,y,as.integer(dims),gamma)
                alpha1  <- pchisq(qchisq(alpha,  p, lower.tail = F)*t(matrix(fit$Sigma2,ncol = q)), p, lower.tail = F)
                pvals   <- pchisq(matrix(rep(fit$Tn,ngamma),nrow = q)*t(matrix(fit$Sigma2,ncol = q)), p, lower.tail = F)
                pvals1  <- pchisq(fit$Tn, p, lower.tail = F)
            }
            else{
                pvals   <- pchisq(fit$Tn*fit$Sigma2^((gamma-1)/gamma), p, lower.tail = F)
            }
        }else{
            pvals <- pchisq(fit$Tn, p, lower.tail = F)
        }
    }else{
        if(is.null(criteria) || criteria=="RBS"){
            if(isbic){
                dims    = c(n,p,q,ngamma)
                fit     = .Call("PVALUES_BIC",x,y,as.integer(dims),gamma)
                alpha1  <- pf(qf(alpha,  df1 = p, df2 = n-p, lower.tail = F)*t(matrix(fit$Sigma2,ncol = q)), df1 = p, df2 = n-p, lower.tail = F)
                pvals   <- pf(matrix(rep(fit$Tn/p, ngamma),nrow = q)*t(matrix(fit$Sigma2, ncol = q)), df1 = p, df2 = n-p, lower.tail = F)
                pvals1  <- pf(fit$Tn/p, p, df1 = p, df2 = n-p, lower.tail = F)
            }
            else{
                pvals   <- pf(fit$Tn*fit$Sigma2^((gamma-1)/gamma)/p, df1 = p, df2 = n-p, lower.tail = F)    
            }
        }else{
            pvals <- pf(fit$Tn/p, df1 = p, df2 = n-p, lower.tail = F)
        }
    }
    
    fit$pvals <- pvals
    indx <- NULL

    if(is.null(criteria)){
        fit$pvfdr <- NA
        fit$signifc <- NA
    }else if(criteria=="BH"){
        sortpv <- sort(pvals, index.return = T)
        thres = alpha*seq(q)/q
        indx = which(sortpv$x<thres)
        if(length(indx)>0){
            j0 = indx[length(indx)]
            fit$pvfdr <- sortpv$x[1:j0]
            fit$signifc <- sortpv$ix[1:j0]
        }else{
            fit$pvfdr <- NA
            fit$signifc <- NA
        }
    }else if(criteria=="RBS"){
        if(isbic){
            indx    <- matrix(rep(pvals1,ngamma),nrow=q) < alpha1
            errhat  <- y - x%*%solve(t(x)%*%x,t(x))%*%y
            bic     <- rep(0,ngamma)
            for(i in 1:ngamma){
                indxi   <- indx[,i]
                tmp     <- y%*%diag(1-indxi) + errhat%*%diag(indxi)
                bic[i]  <- sum(tmp^2)
            }
            bic     <- log(bic/n/q) + log(n*q)*colSums(indx)/n/q
            indx_min = which.min(bic)
            fit$pvfdr <- pvals[,indx_min]
            fit$signifc <- indx[,indx_min]
        }
        else{
            indx <- which(pvals<alpha/q)
            fit$pvfdr <- pvals[indx]
            fit$signifc <- indx
        }
    }else{
        indx <- which(pvals<alpha/q)
        fit$pvfdr <- pvals[indx]
        fit$signifc <- indx
    }

    
    fit
}

