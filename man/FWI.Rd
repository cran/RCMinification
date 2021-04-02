\name{FWI}
\alias{FWI}
\title{Fire Weather Index Series}
\usage{data(FWI)}
\description{
The \code{FWI} list consists of 4 vectors containing daily
Fire Weather Index observations.
}
\format{
  This list contains the following vectors:
  \describe{
    \item{PG2008}{FWI observations from Prince George, BC for 2008}
    \item{PG2009}{FWI observations from Prince George, BC for 2009}
    \item{ED2013}{FWI observations from Edmonton, AB for 2013}
    \item{ED2014}{FWI observations from Edmonton, AB for 2014}
  }
}
\examples{
    RCMTmle(FWI$PG2009[c(100:300)])
}
\keyword{datasets}
