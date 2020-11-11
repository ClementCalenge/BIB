\name{testSpatial}
\alias{testSpatial}
\title{
  Test of the Assumption of Constant Capture Effort Over the Study area.
}
\description{
  
  This function implements the test the assumption of uniform capture
  effort over the core area of the massif.
  
}
\usage{
testSpatial(x, dataSpatial, proba = 0.95, year0 = 2012, nsim = sum(sapply(x, nrow)))
}
\arguments{
  \item{x}{
    An object of class \code{mcmc.list}, containing at least a parameter
    named \code{gammat} corresponding to the proportion of marked
    animals in the population at time t.
  }
  \item{dataSpatial}{
    The data.frame \code{dataSpatial} (see the help page of this dataset).
  }
  \item{proba}{
    Probability used for the credible intervals.
  }
  \item{year0}{
    The year just before the onset of the program.
  }
  \item{nsim}{
    The number of simulations in the test. By default, equal to the
    number of MCMC samples in \code{x}.
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
\value{

  The function returns a data.frame containing, for each combination of
  the years and zones, the number of marked animals with the limits of
  the credible interval on the expected number of marked animals in the
  areas under the model.
  
}
\references{

  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.
  
}
\author{

  Clement Calenge \email{clement.calenge@ofb.gouv.fr}

}
\seealso{
  \code{\link{capturesSpatial}}
}
\examples{

testSpatial(SamplesModelCatalyticYear, capturesSpatial)

}
\keyword{htest}
