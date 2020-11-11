\name{showFOI}
\alias{showFOI}
\title{
  Draws a Boxplot of the Force of Infection
}
\description{
  This function draws a boxplot of the posterior distribution of the
  force of infection by the brucellosis each year between 2013 and
  2018.
}
\usage{
showFOI(x, year0 = 2012)
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
  This function reproduces the figure 4 of Calenge et al. (2020).
}
\value{
  A boxplot.
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
showFOI(SamplesModelCatalyticYear)
}
}
\keyword{hplot}
