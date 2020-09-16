\name{dcorr}
\alias{dcorr}
\title{
  Distance correlation of two multivariates.
}
\description{
  Distance correlation and covariance of two multivariates \eqn{y} and \eqn{x}.
}

\usage{dcorr(y,x)
}

\arguments{
  \item{y}{A \eqn{n\times q} numeric matrix.}
  
  \item{x}{A \eqn{n\times p} numeric matrix.}
}

\value{
  \item{dcor}{The distance correlation, which is an 4-vactor with the dcorr of both \eqn{y} and \eqn{x}, the dcov of \eqn{y}, the dcov of dcorr \eqn{x}, and the dcov of both \eqn{y} and \eqn{x}. dcov denotes the sample distance covariance, and dcorr denotes the sample distance correlation.}
}


\references{
Szekely, G.J. and Rizzo, M.L. (2009). Brownian Distance Covariance, Annals of Applied Statistics, 3(4), 1236-1265.

Szekely, G.J. and Rizzo, M.L. (2009). Rejoinder: Brownian Distance Covariance, Annals of Applied
Statistics, 3(4), 1303-1308.

Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007). Measuring and Testing Dependence by Correlation of Distances, Annals of Statistics, 35(6), 2769-2794.
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

dcor <- dcorr(y,x)

}