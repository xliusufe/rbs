pval <- function(x,y,criteria=NULL,alpha=0.05,gamma=0.15,family="Fdist"){
    n = nrow(x)
    p = ncol(x)
    q = ncol(y)
    if(is.null(p)) p = 1
    if(is.null(q)) q = 1
    
    dims = c(n,p,q)
    fit = .Call("PVALUES",x,y,as.integer(dims))
    if(family=="Chi2"){
        if(is.null(criteria) || criteria=="RBS"){
            pvals   <- pchisq(fit$Tn*fit$Sigma2^(gamma/(1+gamma)), p, lower.tail = F)
        }else{
            pvals <- pchisq(fit$Tn, p, lower.tail = F)
        }
    }else{
        if(is.null(criteria) || criteria=="RBS"){
            pvals   <- pf(fit$Tn*fit$Sigma2^(gamma/(1+gamma))/p, df1 = p, df2 = n-p, lower.tail = F)
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
        indx <- which(pvals<alpha/q)
        fit$pvfdr <- pvals[indx]
        fit$signifc <- indx
    }else{
        indx <- which(pvals<alpha/q)
        fit$pvfdr <- pvals[indx]
        fit$signifc <- indx
    }

    
    fit
}

