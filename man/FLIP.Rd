\name{flip}
\alias{flip}
\title{
  Optimality conditions for the minimization of quadratic 0-1 problems
}
\description{
  Flip procedure for optimality conditions for the minimization of quadratic 0-1 problems, where one-flip, two-flip and hybrid for both are considered. The hybrid flip applies one-flip and two-flip sequentially.
}

\usage{flip(A,b=NULL,x0=NULL,nflip=1)
}

\arguments{
  \item{A}{A \eqn{p}-symmetric matrix.}
  
  \item{b}{A \eqn{p}-vector. Default is zero.}
  
  \item{x0}{The initial value. Default is zero.}
  
  \item{nflip}{An integer \eqn{1, 2, 3} with one-flip if \code{nflip=1}, two-flip if \code{nflip=2}, and hybrid if \code{nflip=3}. Default is \code{nflip=1} corresponding to one-flip.}
}

\value{
  \item{xhat}{The local minimizer.}
  
  \item{obj}{the local minimum.}
}


\references{
Chen, W. (2015). Optimality Conditions for the Minimization of Quadratic 0-1 Problems. SIAM Journal on Optimization, 25(3), 1717-1731.
}

\examples{ 

data(Qd)
Q <- as.matrix(Qd$Q)
fit <- flip(Q,nflip=1)
fit

}