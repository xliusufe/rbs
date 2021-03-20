# RBS
R package "rbs" providing a procedure to select the response variables and estimate regression coefficients simultaneously. It also provides the screening procedure based on the distance correlation, the solution to the quadratic 0-1 programming problems by transferring to k-flipping optimization problems and to continuous quadratic programming problems, and the multiple test procedure including Benjamini-Hochberg and Bonferroni correction.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/rbs")

# Usage

   - [x] [rbs-manual.pdf](https://github.com/xliusufe/rbs/blob/master/inst/rbs-manual.pdf) ---------- Details of the usage of the package.
   - [x] [rbs](https://github.com/xliusufe/RBSSpy) ------------------------ The corresponding `Python` package

# Example
    library(rbs)

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
    y <- cbind(x%*%beta, matrix(0,n,q-q0)) + eps%*%A

    fit <- rbs(x,y,criteria=0)
    fit$delta

    fit <- rbs(x,y)
    fit$delta
    fit$selected

    fit <- rbs_sig(x,y,criteria=0)
    fit$delta


    fit <- rbs_sig(x,y,V,criteria=0)
    fit$delta


    lambda <- seq(0.01, 2, length = 50)
    fit <- rbs_sig(x,y,lambda=lambda)
    fit$delta
    fit$selected

    fit <- rbs_sig(x,y,V,lambda=lambda)
    fit$delta
    fit$selected

    fit <- pval(x,y)
    fit$Tn
    fit$pvals
    fit$pvfdr
    fit$signifc
    
# References

Benjamini, Y. and Hochberg,  Y. (1995). Controlling the False Discovery Rate A Practical and Powerful Approach to Multiple testing. Journal of the Royal Statistical Society: Series B (Methodological). 57(1), 289-300.

Chen, W. and L. Zhang (2010). Global Optimality Conditions for Quadratic 0-1 Optimization Problems. Journal of Global Optimization 46(2), 191-206.

Chen, W. (2015). Optimality Conditions for the Minimization of Quadratic 0-1 Problems. SIAM Journal on Optimization, 25(3), 1717-1731.

Hu, J., Huang, J., Liu, X. and Liu, X. (2020). Response Best-subset Selector for Multivariate Regression. Manuscript.

Li, R., W. Zhong, and L. Zhu (2012). Feature Screening Via Distance Correlation Learning. Journal of the American Statistical Association, 107 (499), 1129-1139.

Szekely, G.J. and Rizzo, M.L. (2009). Brownian Distance Covariance, Annals of Applied Statistics, 3(4), 1236-1265.

Szekely, G.J. and Rizzo, M.L. (2009). Rejoinder: Brownian Distance Covariance, Annals of Applied
Statistics, 3(4), 1303-1308.

Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007). Measuring and Testing Dependence by Correlation of Distances, Annals of Statistics, 35(6), 2769-2794.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
