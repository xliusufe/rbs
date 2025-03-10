\name{sirs}
\alias{sirs}
\title{
  Sure independence ranking and screening (SIRS)
}
\description{
  A model-free screening method, proposed by Zhu et al. (2011), 
  which combines a soft cutoff value and a hard cutoff.
}

\usage{sirs(Y, X, standardize_X = TRUE, N = NULL, d = 10, ntop = 10)
}

\arguments{
  \item{Y}{An \eqn{n}-dimensional numeric vector.}
  
  \item{X}{An \eqn{n\times p} numeric matrix.}
  
  \item{standardize_X}{Logical flag for \eqn{X} standardization, prior to performing screening.
  Default is \code{standardize_X=TRUE}.}

  \item{N}{An integer specified by the user for the hard threshold rule, representing the number of 
  top N ranked variables to retain after sorting. The default value is \code{NULL}, 
  which corresponds to \eqn{\lfloor n / \log(n) \rfloor} (floor function).}

  \item{d}{An integer specified by the user for the soft threshold rule, representing the 
  dimension of the auxiliary variables. The default value is 10.}

  \item{ntop}{An integer, which is integer that the indices of the top \code{ntop} most correlated variables will be output.}
}

\value{
  \item{indn}{The indices of the top \code{ntop} most correlated variables. 
  If \code{ntop} exceeds the number of variables selected by the combined soft and hard threshold rules, 
  the indices of all variables selected by the threshold rules are returned. 
  Otherwise, the indices of the top \code{ntop} variables are returned.}
}

\references{
  Zhu, L. P., Li, L., Li, R. and Zhu, L. X. (2011). Model-free feature screening for ultrahigh-dimensional data. 
  Journal of the American Statistical Association, 
  106(496), 1464-1475.
}

\examples{ 
n   <- 100
p   <- 200

beta <- c(rep(0.5, 20), rep(0, p - 20))
eps <- rnorm(n, mean = 0, sd = 1)

X <- matrix(rnorm(n*p),n,p)
Y <- X \%*\% beta + eps

fit <- sirs(Y,X)
}