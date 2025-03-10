\name{pcorr}
\alias{pcorr}
\title{
  Pearson correlation coefficient
}
\description{
  Pearson correlation coefficient vector between an n-dimensional vector \eqn{Y} and an n-by-p matrix \eqn{X}.
}

\usage{
  pcorr(Y,X)
}

\arguments{
  \item{Y}{An \eqn{n}-dimensional numeric vector.}
  
  \item{X}{An \eqn{n\times p} numeric matrix.}
}

\value{
  \item{pcor}{A \eqn{p}-dimensional vector, where the i-th element corresponds to the Pearson correlation 
  coefficient between \eqn{Y} and the i-th column of \eqn{X}.}
}

\examples{ 
n   <- 100
p   <- 200

beta <- c(rep(0.5, 20), rep(0, p - 20))
eps <- rnorm(n, mean = 0, sd = 1)

X <- matrix(rnorm(n*p),n,p)
Y <- X \%*\% beta + eps

pcor <- pcorr(Y,X)

}