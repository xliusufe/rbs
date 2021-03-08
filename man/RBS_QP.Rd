\name{rbs_qp}
\alias{rbs_qp}
\title{
  RBSS with considering covariance of responses based on continuous quadratic programming.
}
\description{
  Select the response variables and estimate regression coefficients simultaneously for multivariate linear regression with considering covariance of responses, in which the quadratic 0-1 programming problems are transferred to continuous quadratic programming problems.
}

\usage{rbs_qp(x,y,V=NULL,gamma=1.5,lambda=NULL,criteria=2,tau=1)
}

\arguments{
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{y}{A \eqn{n\times q} response matrix.}
  
  \item{V}{A user-specified \eqn{q\times q} precision matrix. A estimator is provided if \code{V=NULL}. Default is \code{V=NULL}.}
  
  \item{gamma}{A positive separating parameter \eqn{\gamma}. Default is \code{gamma=1.5}.}
  
  \item{lambda}{A user-specified sequence of \eqn{\lambda} values.  By default,
    a sequence of values of length \code{nlambda} is computed, equally
    spaced on the scale.}
    
  \item{criteria}{The criteria to be applied to select parameters.  Either \code{AIC} if \code{criteria=1}, \code{BIC} (the default) if \code{criteria=2}, 
                or \code{GCV} if \code{criteria=3}. There is no selection if \code{criteria=0}, in which case \code{lambda} should be a number.}   
                
  \item{tau}{A constant to adjust \code{AIC} creteria. Default is \code{tau=1}.}
}

\value{
  \item{delta}{The estimation of the \eqn{\delta}.}
  \item{theta}{The estimation of the \eqn{\theta}.}
  \item{rss}{Residual sum of squares (RSS) without the selection of tuning parameters.}
  
  \item{deltapath}{The estimation of the \eqn{\delta} with the selection of tuning parameters.}
  \item{bic}{The \code{AIC} or \code{BIC} or \code{GCV} with the selection of tuning parameters.}
  \item{selected}{The index of \eqn{\lambda} corresponding to \code{lambda_opt} with the selection of tuning parameters.}
}


\references{
Chen, W. and L. Zhang (2010). Global Optimality Conditions for Quadratic 0-1 Optimization Problems. Journal of Global Optimization 46(2), 191-206.

Hu, J., Huang, J., Liu, X. and Qiu F. (2020). Response Best-subset Selector for Multivariate Regression. Manuscript.
}

\examples{ 
n   <- 200
p   <- 5
q   <- 10
q0  <- 5

Sigma <- matrix(0,q,q)
for(i in 1:q) for(j in 1:q) Sigma[i,j]=0.5^(abs(i-j))
A <- chol(Sigma)
V <- solve(Sigma)

beta <- matrix(runif(p*q0),p,q0)
eps <- matrix(rnorm(n*q),n,q)

x <- matrix(rnorm(n*p),n,p)
y <- cbind(x\%*\%beta, matrix(0,n,q-q0)) + eps\%*\%A

fit <- rbs_sig(x,y,lambda=0.4)
fit$delta

fit <- rbs_sig(x,y,V,lambda=0.4)
fit$delta

lambda <- seq(0.01, 2, length = 50)
fit <- rbs_sig(x,y,lambda=lambda)
fit$delta
fit$selected

fit <- rbs_sig(x,y,V,lambda=lambda)
fit$delta
fit$selected
}