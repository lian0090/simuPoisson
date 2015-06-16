\name{map}
\alias{map}
\docType{data}
\title{
map positions for markers
}
\description{
map positions for makrers with marker name, chromsome number, position on individual chromosomes and cululated position on all chromosomes
}
\usage{data(map)}
\format{
  A data frame with 2911 observations on the following 4 variables.
  \describe{
  	\item{\code{mkname}}{a factor with all the marker names}
    \item{\code{chr}}{a numeric vector}
    \item{\code{pos}}{a numeric vector of map positions on each individual chromosome in centimorgan. Should all start with 0}
    \item{\code{cumpos}}{Optional: a numeric vector of cumulated map positions on all chromosomes }
  }
}

\source{
Monsanto Breeding data
}
\examples{
data(map)
str(map) 
}
