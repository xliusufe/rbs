\name{rbsrp}
\alias{rbsrp}
\title{
  RBS without covariance of responses based on the random projection approach
}
\description{
  Select the response variables and estimate regression coefficients simultaneously for multivariate linear regression without covariance of responses based on the random projection approach.
}

\usage{rbsrp(x,y,rho=0.4,Pk=NULL,gamma=1.5, lambda=NULL,criteria=2,tau=1)
}

\arguments{
  \item{x}{A \eqn{n\times p} Numeric design matrix for the model.}
  
  \item{y}{A \eqn{n\times q} Response matrix.}

  \item{rho}{The given ratio \eqn{\rho=\lim \frac{k}{n}}, a scale in the unit inteval. Default is \code{rho = 0.4}.}

  \item{Pk}{The preset random matrix. Default is \code{Pk = NULL}, in the case of which each entry is drawn independently from the standard normal distribution.}
  
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
  
  \item{deltapath}{The estimation path of the \eqn{\delta} with the selection of tuning parameters.}
  \item{bic}{The \code{AIC} or \code{BIC} or \code{GCV} with the selection of tuning parameters.}
  \item{selected}{The index of \eqn{\lambda} corresponding to \code{lambda_opt} with the selection of tuning parameters.}
  \item{rho}{The input ratio \eqn{\rho=\lim \frac{k}{n}}.}
}


\references{
Hu, J., Huang, J., Liu, X. and Liu, X (2020). Response Best-subset Selector for Multivariate Regression. Manuscript.
}

\examples{ 
n   <- 200
p   <- 5
q   <- 10
q0  <- 5

beta <- matrix(runif(p*q0),p,q0)
eps <- matrix(rnorm(n*q),n,q)

x <- matrix(rnorm(n*p),n,p)
y <- cbind(x\%*\%beta, matrix(0,n,q-q0)) + eps

fit <- rbsrp(x,y,lambda=0.4)
fit$delta

lambda <- seq(0.01, 2, length = 50)
fit <- rbsrp(x,y,lambda=lambda)
fit$delta
fit$selected

}