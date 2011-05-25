\name{MEWMA-package}
\alias{MEWMA-package}
\alias{MEWMA}
\docType{package}
\title{Multivariate Exponentially Weighted Moving Average (MEWMA) Control Chart}

\description{
MEWMA package perform the MEWMA Control Chart according to Lowry et al.(1992).
Allow to monitor simultaneously multiple quality characteristics.
MEWMA is an extension of the univariate EWMA chart. 
}
\details{
\tabular{ll}{
Package: \tab mewma\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-25\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
The mewma function have the following structure, mewma(x,n,lambda); where x must be a matrix, n(subgroup size) 
an integer and lambda(smoothing constant) a value between 0.1 and 0.9.   
The cases of individual observations the covariance matrix is estimated according to Holmes and Mergen(1993).
The package use the UCL proposed by Bodden and Rigdon(1999) and only for the followings parameters:
quality characteristics between 2 and 10, smoothing constant(lambda)=0.1,0.2,0.3,...,0.9 and ARL=200.

}
\author{
Edgar Santos Fernandez

Maintainer: Edgar Santos Fernandez <edgar.santos@etecsa.cu>

}

\references{

Bodden,Kevin M. and Rigdon,Steven E. (1999)
"A Program for Approximating the In-Control ARL for the MEWMA Chart",Journal of Quality Technology
31.

Holmes, Donald and Mergen, A. Erhan (1993)
"Improving the Performance of the T2 Control Chart",Quality Engineering
5.

Lowry, Cynthia A. and Woodall, William H. and Champ, Charles W. and Rigdon, Steven E.(1992)
"A Multivariate Exponentially Weighted Moving Average Control Chart",Journal of Quality Technology
34.

 Montgomery,Douglas C.(2005)
"Introduction to Statistical Quality Control",John Wiley & Sons
5 edition.

}

\keyword{ package }
\seealso{
For more details see Lowry et al.(1992)
}
\examples{
# constructing a MEWMA chart with subgroup size(n)=2 and lambda=0.1
x<-rgamma(300,1,3); x<-matrix(x,nrow=60)
mewma(x,2,0.1)

# constructing a MEWMA chart for individual observations
x<-rgamma(500,1,1); x<-matrix(x,nrow=100)
mewma(x,1,0.2)

# constructing a MEWMA chart and saving the t2 statistics
x<-rbeta(400,1,2); x<-matrix(x,nrow=100)
a<-mewma(x,5,0.1)$t2;print(a)

# constructing a MEWMA chart and saving the covariance matrix
x<-rbeta(400,1,2); x<-matrix(x,nrow=100)
b<-mewma(x,5,0.1)$covariance;print(b)
}
