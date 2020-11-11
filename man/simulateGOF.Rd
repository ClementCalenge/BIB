\name{simulateGOF}
\alias{simulateGOF}
\alias{print.gofBIB}
\title{
  Implements the Goodness of Fit Tests of the Models
}
\description{
  Implements the Goodness of Fit Tests of the Models
}
\usage{
simulateGOF(x, dataset, nsim = sum(sapply(x, nrow)), proba = 0.8)

\method{print}{gofBIB}(x, \dots)

}
\arguments{
  \item{x}{
    For \code{simulateGOF}, an object of class \code{mcmc.list},
    corresponding to the model \code{ModelCatalyticYear} (see the help
    page of this model). Can only be used currently with the object
    \code{SamplesModelCatalyticYear} For \code{print.gofBIB} an object
    of class \code{gofBIB}. 
  }
  \item{dataset}{
    The dataset used to fit the model.
  }
  \item{nsim}{
    The number of simulations used to fit the model.
  }
  \item{proba}{
    The probability of the credible intervals used to assess the fit.
  }
  \item{\dots}{
    additional arguments that can be passed to the function \code{print}
  }
}
\details{

  This function implements the approach recommended by Gelman and Meng
  (1996). Each MCMC iteration generated a sampled value the vector of
  parameters. For each simulated value of this vector, this function
  simulates a replication of the captures, bacterial cultures and
  censuses datasets. This function then compares summary statistics
  calculated on the observed datasets to the distribution of these
  summary statistics calculated on the simulated datasets. More
  specically, this function focuses on the following statistics: (i) the
  number of marked females detected during each day of census in the
  censuses dataset, (ii) the proportion of unmarked females detected as
  seropositive each year in the seroprevalence dataset, (iii) the same
  proportion calculated only on marked females, (iv) the log-titer of
  the CFT for each seropositive female.

}
\value{

  An object of class \code{gofBIB}.
  
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
}
\examples{

\dontrun{

simulateGOF(SamplesModelCatalyticYear, BruData)

}
}
\keyword{model}
