\name{robustnessPitWit}
\alias{robustnessPitWit}
\title{
  Tests the Robustness of the ModelCatalyticYear for Various Priors
}
\description{
  This function fits the model ModelCatalyticYear with different beta priors
  for the prevalence among unmarked animals (\code{pitu}) and for the
  proportion of seropositive animals that are actively infected
  (\code{wit}).
}
\usage{
robustnessPitWit(be1p = 1, be2p = 1, be1w = 1, be2w = 1, dataset, niter = 5000, thin = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{be1p}{
    First shape parameter for \code{pitu}
  }
  \item{be2p}{
    Second shape parameter for \code{pitu}
  }
  \item{be1w}{
    First shape parameter for \code{wit}
  }
  \item{be2w}{
    Second shape parameter for \code{wit}
  }
  \item{dataset}{
    Dataset required to fit the model (see \code{?BruData})
  }
  \item{niter}{
    Number of MCMC iterations to fit the model
  }
  \item{thin}{
    Thinning interval in the MCMC
  }
}
\value{
  An object of class \code{mcmc.list} containing the MCMC samples of the model.
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
  \code{\link{BruData}} for a description of the structure required for
  the parameter \code{dataset}, \code{\link{ModelCatalyticYear}} for a
  description of the model.
}
\examples{

\dontrun{

## Simulate the model with a beta(2,2) prior for pitu and beta(5,5)
## prior for wt
robustnessPitWit(2, 2, 5, 5, BruData)

}


}
\keyword{htest }
