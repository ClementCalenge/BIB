\name{getTable2and3}
\alias{getTable2and3}
\title{
  Main Parameters of the Model Predicting the Seroprevalence of
  Brucellosis in the Ibex Population.
}
\description{
  This function returns the table 2 and 3 of the paper containing the
  posterior mode and credible intervals on the main parameters of the
  model.
}
\usage{
getTable2and3(x, years = 2013:2018, proba = 0.9, table2 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{mcmc.list}, containing at least a parameter
    named \code{pitu}, \code{pitm}, \code{pit} and \code{gammat}
    corresponding to the proportion of marked animals in the population
    at time t.
}
  \item{years}{
    The years of interest.
}
  \item{proba}{
    The probability of the credible intervals.
}
  \item{table2}{
    Logical. Whether the function will return the posterior mode and
    credible intervals of table 2 (if \code{TRUE}), or of table 3 (if
    \code{FALSE}). 
}
}
\details{

  The table 2 reports results of the overall seroprevalence,
  seroprevalence of marked animals, seroprevalence of unmarked animals,
  proportion of marked animals.

  The table 3 reports results of the proportion of seropositive animals
  with active infection, the overall proportion of animals with active
  infection.
  
}
\value{

  A data.frame.
  
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
  \code{\link{SamplesModelCatalyticYear}} for a description of the model
  fit to the data.
}
\examples{

## Table 2
getTable2and3(SamplesModelCatalyticYear)

## Table 3
getTable2and3(SamplesModelCatalyticYear, table2 = FALSE)

}
\keyword{models}
