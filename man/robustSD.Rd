\name{robustSD}
\alias{robustSD}
\title{
Tatum's Robust Estimate of the Standard Deviation
}
\description{
Standard deviation estimate which is insensitive to outliers and
random trends.  
}
\usage{
robustSD(x)
}
\arguments{
  \item{x}{A numeric vector.} 
}
\value{
a numeric vector of length 1. 
}
\author{
L. Han
}
\references{Tatum, L.G. (1997) Robust Estimation of the Process Standard
Deviation for Control Charts. Journal of the American Statistical Association
39, pp 127-141.  }

\examples{
robustSD(EuStockMarkets[,1])
}
\keyword{statistics}
