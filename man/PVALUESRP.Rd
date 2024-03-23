\name{pvalrp}
\alias{pvalrp}
\title{
  P-values for F-test of the separate responses based on the random projection approach
}
\description{
   P-values for F-test of the separate responses for the multivariate linear regression models based on the random projection approach.
}

\usage{pvalrp(x,y,rho=0.4,Pk=NULL,criteria=NULL,alpha=0.05,gamma=1.15,
                                              family="Fdist",isbic=FALSE)
}

\arguments{
  \item{x}{A \eqn{n\times p} Numeric design matrix for the model.}
  
  \item{y}{A \eqn{n\times q} Response matrix.}

  \item{rho}{The given ratio \eqn{\rho=\lim \frac{k}{n}}, a scale in the unit inteval. Default is \code{rho = 0.4}.}

  \item{Pk}{The preset random matrix. Default is \code{Pk = NULL}, in the case of which each entry is drawn independently from the standard normal distribution.} 

  \item{criteria}{A criteria to select important variables by a significant level. No corrections if \code{criteria=NULL}, RBS
procedure if \code{criteria="RBS"}, Benjamini-Hochberg
procedure if \code{criteria="BH"}, and Bonferroni correction if \code{criteria="Bonf"}.}
  
  \item{alpha}{A prespecified level.}
  
  \item{gamma}{A positive separating parameter \eqn{\gamma} if RBS procedure is used. Default is \code{gamma=1.15}.}
  
  \item{family}{A string representing one of the built-in families, by which P-values are calculated. F-test is used if \code{family="Fdist"} with the first degrees of freedom \eqn{p} and the second degrees of freedom \eqn{n-p}, and \eqn{\chi^2}-test is used if \code{family="Chi2"} with degrees of freedom \eqn{p}. Default is \code{family="Fdist"} (F-test).}
   
   \item{isbic}{A logical flag. The \code{BIC} criteria is used (\code{TRUE}) or not (\code{default = FALSE}).}
}

\value{
  \item{Tn}{Values of test statistics.}
  \item{Sigma2}{Estimator of the marginal response variance.}
  \item{pvals}{P-values.}
  \item{pvfdr}{The P-values corresponding to selected variables.}
  \item{signifc}{The indices corresponding to selected variables.}
  \item{rho}{The input ratio \eqn{\rho=\lim \frac{k}{n}}.}
}


\references{
Benjamini, Y. and Hochberg,  Y. (1995). Controlling the False Discovery Rate A Practical and Powerful Approach to Multiple testing. Journal of the Royal Statistical Society: Series B (Methodological). 57(1), 289-300.

Hu, J., Huang, J., Liu, X. and Liu, X. (2020). Response Best-subset Selector for Multivariate Regression. Manuscript.
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

fit <- pvalrp(x,y)

fit$Tn
fit$pvals
fit$pvfdr
fit$signifc

}