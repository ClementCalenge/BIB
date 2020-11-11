\name{BIB-package}
\alias{BIB-package}
\alias{BIB}
\docType{package}
\title{
\packageTitle{BIB}
}
\description{
\packageDescription{BIB}
}
\details{
  
  This package is the companion package to the paper of Calenge et
  al. (in prep.) describing the estimation of the prevalence of the
  brucellosis in Alpine ibex in the Bargy (France).  It contains all the
  data and functions to reproduce the modelling approach presented in
  this paper.  A vignette details the whole process and model, but we
  summarize it here.
  
  The dataset \code{BruData} is the dataset used to fit the model.  The
  model is fit with the functions of the package \code{rjags}.  The JAGS
  model is coded in the object \code{ModelCatalyticYear}.  Other
  alternative models studied in the appendices of the paper are stored
  in the objects \code{ModelCatalytic} (catalytic model),
  \code{ModelConstPrev} (constant prevalence), and \code{ModelDensity}
  (force of infection proportional to the density).  The result of the
  fit of these models is stored in the objects
  \code{SamplesModelCatalyticYear}, \code{SamplesModelCatalytic},
  \code{SamplesModelConstPrev}, and \code{SamplesModelDensity}.
  
  The package provided several functions that can be used to generate
  the tables and figures used in the paper: \code{getTable2and3}
  (generates the table 2 and 3 of the paper), \code{showFOI} (show the
  changes in force of infection, as in Fig. 4 of the paper),
  \code{showSummary} (summarizes the results, as in Fig. 3 of the
  paper), \code{simulateGOF} (tests the goodness of fit of the model),
  \code{totalPrevalence} (generates the dataset that can be used to
  compare alternative models).

  Note that the test of the assumption of constant capture
  effort over the core area of the massif can be carried out with the
  function \code{testSpatial} carried out on the dataset
  \code{capturesSpatial}. 

}
\author{
\packageAuthor{BIB}

Maintainer: \packageMaintainer{BIB}
}
\references{
  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.
}
\keyword{ package }
\examples{

## All relevant examples are presented on the help pages
## of the functions described in the section details

}
