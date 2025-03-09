\name{pcorr}
\alias{pcorr}
\title{
  Pearson correlation coefficient
}
\description{
  Pearson correlation coefficient vector between an n-dimensional vector \eqn{y} and an n-by-p matrix \eqn{x}.
}

\usage{
  pcorr(y,x)
}

\arguments{
  \item{y}{An \eqn{n}-dimensional numeric vector.}
  
  \item{x}{An \eqn{n\times p} numeric matrix.}
}

\value{
  \item{pcor}{A \eqn{p}-dimensional vector, where the i-th element corresponds to the Pearson correlation 
  coefficient between \eqn{y} and the i-th column of \eqn{x}.}
}

\examples{ 
n   <- 100
p   <- 200

beta <- c(rep(0.5, 20), rep(0, p - 20))
eps <- rnorm(n, mean = 0, sd = 1)

x <- matrix(rnorm(n*p),n,p)
y <- x %*% beta + eps

pcor <- pcorr(y,x)

}