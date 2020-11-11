\name{totalPrevalence}
\alias{totalPrevalence}
\title{
  Calculates the Overall Prevalence of the Brucellosis for Each Year of
  the Study Period
}
\description{
  
  This function estimates the point estimate and the credible intervals
  of the overall prevalence of the Brucellosis for each year of the
  study period.

}
\usage{
totalPrevalence(x, limits = c(0.1, 0.9), years = 2013:2018)
}
\arguments{
  \item{x}{
    An object of class \code{mcmc.list},
    corresponding to the model \code{ModelCatalyticYear} (see the help
    page of this model).  Can only be used currently with the object
    \code{SamplesModelCatalyticYear}. 
  }
  \item{limits}{
    The limits of the credible intervals.
  }
  \item{years}{
    The years of the monitoring.
  }
}
\value{
  
  A data.frame giving the point estimate and limits of the credible
  interval on the overall prevalence.

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
    \code{\link{SamplesModelCatalyticYear}}, \code{\link{ModelCatalyticYear}}
}
\examples{

totalPrevalence(SamplesModelCatalyticYear)

}
\keyword{model}
