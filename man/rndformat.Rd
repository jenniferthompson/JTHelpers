% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rndformat.R
\name{rndformat}
\alias{rndformat}
\title{Round and format a number to the same number of digits.}
\usage{
rndformat(x, digits = 2)
}
\arguments{
\item{x}{Numeric; value to be rounded and formatted}

\item{digits}{Numeric; number of digits to round and format x to}
}
\value{
Character string of x, rounded and formatted to the same number of digits.
}
\description{
Round and format a number to the same number of digits.
}
\examples{
rndformat(1.9727, digits = 3)
rndformat(1.2, digits = 3)

}

