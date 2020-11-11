\name{capturesSpatial}
\alias{capturesSpatial}
\docType{data}
\title{
  Dataset Used to Test the Assumption of Constant Capture Effort
}
\description{
  This data.frame contains the data used to test the assumption of
  uniform capture effort over the core area of the massif.
}
\usage{data("capturesSpatial")}
\format{
  This data.frame contains, for each combination of zone and year, the
  following variables:
  \describe{
    \item{\code{Year}}{a numeric vector containing the year}
    \item{\code{Ntot}}{a numeric vector containing the number of animals
      seen during the censuses.}
    \item{\code{Nmarked}}{a numeric vector containing the number of
      marked animals seen during the censuses.}
    \item{\code{Zone}}{a factor with levels \code{Southwest} \code{Northeast}}
  }
}
\details{

  We divided our study area into two halves: a northeastern part and a
  southwestern part. We used our model to test whether the proportion of
  marked animals observed in these two areas could have been expected
  under our non-spatialized model.  This test is implemented by the
  function \code{testSpatial}, and the examples of this help page
  describe how to use it.

}
\source{
  The original Dataset has been kindly provided by XXXXXXXXX
}
\references{

  
  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.

}
\examples{

## see ?testSpatial for examples

}
\keyword{datasets}
