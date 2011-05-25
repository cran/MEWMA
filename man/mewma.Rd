\name{mewma}
\alias{mewma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Exponentially Weighted Moving Average (MEWMA) Control Chart
}

\description{
mewma function perform the MEWMA Control Chart according to Lowry et al.(1992).
Allow to monitor simultaneously multiple quality characteristics.
MEWMA is an extension of the univariate EWMA chart. 
}
\usage{
mewma(x, n, lambda)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x must be a matrix
}
  \item{n}{
n(subgroup size) an integer 
}
  \item{lambda}{
smoothing constant(lambda) must take a value of the followings (0.1,0.2,0.3,...,0.9)
}
}
\details{
For more details see Lowry et al.(1992)
}
\value{

\item{x}{Matrix of the quality characteristics}
\item{n}{Subgroup size}
\item{lambda }{Smoothing constant}
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

}
\author{
Edgar Santos Fernandez
}
\note{
If the values of n and lambda are different of the specified in arguments section, an error message is displayed.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
