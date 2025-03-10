\name{sis}
\alias{sis}
\title{
  Sure independence screening (SIS)
}
\description{
  Screening procedure based on Pearson correlation coefficient of \eqn{Y} and \eqn{X}.
}

\usage{sis(Y, X, ntop=10)
}

\arguments{
  \item{Y}{An \eqn{n}-dimensional numeric vector.}
  
  \item{X}{An \eqn{n\times p} numeric matrix.}
  
  \item{ntop}{An integer specifying the number of top correlated variables (\code{ntop}) whose indices will be output.}
}

\value{
  \item{pcor}{The whole Pearson correlation.}
  
  \item{indn}{The indices of the top \code{ntop} most correlated variables.}
}


\references{
  Fan, J. and Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space. 
  Journal of the Royal Statistical Society Series B: Statistical Methodology, 
  70(5), 849-911.
}

\examples{ 
n   <- 100
p   <- 200

beta <- c(rep(0.5, 20), rep(0, p - 20))
eps <- rnorm(n, mean = 0, sd = 1)

X <- matrix(rnorm(n*p),n,p)
Y <- X \%*\% beta + eps

fit <- sis(Y,X)
}