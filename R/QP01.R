
biQP <- function(Q, b, e){
  d = dim(Q)[1]
  qs = diag(Q)
  Q_u <- Q_ <- Q
  diag(Q_) = 0
  
  # adjust
  u = ((1 - e) * max(apply(abs(Q_), 1, sum)) - min(eigen(Q_)$values)) /
    e + 10 ^ (-8)
  
  diag(Q_u) = u
  b_u = 0.5 * qs + b - 0.5 * u
  mat = c(rep(c(1, -1, rep(0, 2 * d)), d - 1), 1, -1)
  Amat <- t(matrix(mat, 2 * d))
  bvec = rep(c(0, -1), d)
  solve.QP(Q_u, -b_u, Amat, bvec)$solution
}

ridge <- function(y, x, lambda2) {
  n = dim(x)[1]
  p = dim(x)[2]
  xty = t(x) %*% y
  if (n >= p) {
    U = chol(t(x) %*% x / n + 2 * lambda2 * diag(1, p))
    L = t(U)
    bhat = solve(U, solve(L, xty / n))
  } else{
    rho = 2 * n * lambda2
    U = chol(diag(1, n) + x %*% t(x) / rho)
    L = t(U)
    bhat = xty / rho - t(x) %*% solve(U, solve(L, x %*% xty)) / rho ^ 2
  }
  bhat
}

Alg<-function(x,y,v=NULL,ga,r){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  
  A <- B <- matrix(0,q,q)
  
  px=x%*%solve(t(x)%*%x)%*%t(x)
  if(is.null(v))  v = solve( (t(y)%*%(diag(n)-px)%*%y)/(n-p) )
  
  for(i in 1:q){
    for(j in i:q){
      A[i,j] = v[i,j]*(y[,i]%*%(diag(n)-px)%*%y[,j])
      B[i,j] = v[i,j]*y[,i]%*%px%*%y[,j]
      if(i<j){
        A[j,i] = A[i,j]
        B[j,i] = B[i,j]
      }
    }
  }
  diag(B)=(diag(B))^ga
  cf=biQP(A+r*B,-r*B%*%rep(1,q),0.9)
  ifelse(cf<0.5,0,1)
}

rbs_qp <- function(x,y,V=NULL,gamma=1.5,lambda=NULL,criteria=2,tau=1){
  n = nrow(y)
  q = ncol(y)
  p = ncol(x)
  if(is.null(p)) p = 1
  b0 = ridge(y,x,0)
  if(criteria==0){
    dmat <- Alg(x,y,V,gamma,lambda)
    fit <- list()
    fit$delta <- dmat
    fit$theta <- ridge(y%*%diag(as.vector(dmat)),x,0)
    fit$rss <- base::norm(y-x%*%fit$theta%*%diag(as.vector(dmat)),"F")
  }
  else{
    s <- length(lambda)
    dmat = matrix(0,s,q)
    bic = rep(0,s)
    for(i in 1:s){
      dmat[i,] <- tmp <- Alg(x,y,V,gamma,lambda[i])
      loss = base::norm(y-x%*%b0%*%diag(tmp),"F")
      df = sum(tmp!=0)
      if(criteria==1)
        bic[i] = log(loss/n/q) + 2*((p+1)*df)^tau/n/q
      else if(criteria==2)
        bic[i] = log(loss/n/q) + log(n*q)*(p+1)*df/n/q
      else
        bic[i] = loss*n*q /(n*q-df)/(n*q-df)
    }
    id = which.min(bic)
    
    fit <- list()
    fit$delta <- dmat[id,]
    fit$theta <- ridge(y%*%diag(as.vector(dmat)),x,0)
    fit$deltapath <- dmat
    fit$bic <- bic
    fit$selected <- id
  }
  return(fit)
}

