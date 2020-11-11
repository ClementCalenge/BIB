\name{ModelCatalytic}
\alias{ModelCatalytic}
\alias{ModelCatalyticYear}
\alias{ModelConstPrev}
\alias{ModelDensity}
\docType{data}
\title{
  JAGS Models Compared in the Paper
}
\description{
  These datasets are character strings containing the 
}
\usage{

data("ModelCatalytic")

data("ModelCatalyticYear")

data("ModelConstPrev")

data("ModelDensity")

}
\format{
  The object is a character string containing the JAGS code for the models
}
\details{
  
  These models correspond to the four models compared in the paper.  The
  main model discussed in the paper corresponds to the model
  \code{ModelCatalyticYear}.  The model \code{ModelCatalytic}
  corresponds to the catalytic infection model.  The model
  \code{ModelConstPrev} corresponds to the assumption of constant
  prevalence among recaptured animals.  The model \code{ModelDensity}
  corresponds to a force of infection dependent on the density of the
  species.
  
}
\references{

  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.

}
\examples{

cat(ModelCatalyticYear)


}
\keyword{datasets}
