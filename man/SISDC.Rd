\name{sisdc}
\alias{sisdc}
\title{
  Screening procedure based on the distance correlation.
}
\description{
  Screening procedure based on the distance correlation of two multivariates \eqn{y} and \eqn{x}.
}

\usage{sisdc(y, x, d=1, ntop=10)
}

\arguments{
  \item{y}{A \eqn{n\times q} numeric matrix.}
  
  \item{x}{A \eqn{n\times p} numeric matrix.}
  
  \item{d}{An integer. Screening variable \eqn{y} if \code{d=1}, and Screening variable \eqn{x} if \code{d=2}.}
  
  \item{ntop}{An integer, which is integer that the indices of the top \code{ntop} most correlated variables will be output.}
}

\value{
  \item{dcor}{The whole distance correlation.}
  
  \item{indn}{The indices of the top \code{ntop} most correlated variables.}
}


\references{
Li, R., W. Zhong, and L. Zhu (2012). Feature Screening Via Distance Correlation Learning. Journal of the American Statistical Association, 107 (499), 1129-1139.

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

fit <- sisdc(y,x)
fit
}