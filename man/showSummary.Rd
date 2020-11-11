\name{showSummary}
\alias{showSummary}
\title{
  Summary Plot Showing the Changes with Time of Key Parameters of the
  Brucellosis Infection Process in the Bargy Mountains.
}
\description{
  This function shows the changes with time of the overall prevalence of
  the disease, the proportion of marked animals in the population, the
  proportion of seropositive animals with an active infection, and the
  proportion of all animals with an active infection.
}
\usage{
showSummary(x, year0 = 2012)
}
\arguments{
  \item{x}{
    An object of class \code{mcmc.list},
    corresponding to the model \code{ModelCatalyticYear} (see the help
    page of this model).  Can only be used currently with the object
    \code{SamplesModelCatalyticYear}. 
  }
  \item{year0}{
    The year before the onset of the study.
  }
}
\details{
  This function reproduces the figure 3 of Calenge et al. in prep.
}
\value{
  A multipanel plot.
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
\dontrun{
showSummary(SamplesModelCatalyticYear)
}
}
\keyword{hplot}
