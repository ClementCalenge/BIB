\name{overlapPriorPost}
\alias{overlapPriorPost}
\title{
  Calculate the Overlap between the Prior and the Posterior
}
\description{
  This function implements the diagnostic of Garrett and Zeger (2000) to
  assess the overlap between the prior and the posterior.
}
\usage{
overlapPriorPost(x, param, prior = "dnorm(z,0,sqrt(1000))",
                 from = -10, to = 10, n = 10000)
}
\arguments{
  \item{x}{
    An object of class \code{mcmc.list}
  }
  \item{param}{
    a character string containing the name of a parameter in \code{x}.
  }
  \item{prior}{
    a character string indicating the form of the prior for the
    parameter.  Note that the parameter should be named \code{z} in this
    expression.
  }
  \item{from,to,n}{
    the values of \code{param} over which the prior and posterior should
    be compared (passed to the function \code{density}.
  }
}
\value{
  A numeric value corresponding to the diagnostic.
}
\references{
  
  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.
  
  Garrett, E. and Zeger, S. 2000. Latent class model
  diagnosis. Biometrics, 56, 1055-1067.

}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}

\seealso{
  \code{\link{SamplesModelCatalyticYear}}
}
\examples{

## example with the parameter beta0
overlapPriorPost(SamplesModelCatalyticYear, "beta0")

}
\keyword{htest}
