\name{Globaltemps}
\alias{Globaltemps}
\title{Global Average Temperature Changes}
\usage{data(Globaltemps)}
\description{
Global average temperatures are recorded in terms of number of Celsius degrees above a 
baseline temperature from 1880 to 2016.  The baseline temperature is 
the average temperature for the year 1990.
}
\format{
  A numeric vector 
}
\examples{
    temps <- ts(Globaltemps, start = 1880, end = 2016) 
    ts.plot(temps, ylab =  "Change in Temperature")
}
\keyword{datasets}
