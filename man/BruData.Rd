\name{BruData}
\alias{BruData}
\docType{data}
\title{
  Simulated Dataset with Properties Similar to the Dataset  Used to Fit the Model in Calenge et al. (in prep.)
}
\description{
  This dataset has been simulated to replace the sensitive dataset used
  by Calenge et al. (in prep). It contains the captures dataset,
  censuses dataset, and bacterial cultures dataset used to fit the model
  of Calenge et al. (in prep.).
}
\usage{data("BruData")}
\format{
  A list containing the following variables:
  \describe{
    \item{CensusesNmark}{Number of marked animals during each census day.}
    \item{CensusesNind}{Total number of marked animals during each census day.}
    \item{CensusesYear}{Year corresponding to each census day.}
    \item{CensusesNdays}{Number of census days.}
    \item{CapturesN}{Number of captured animals}
    \item{CapturesYear}{Year corresponding to each captured animal}
    \item{CapturesSeroplus}{Whether each captured animal is seropositive.}
    \item{CaptureNumber}{Number of captured animals during each year.}
    \item{RecapturesN}{Number of recaptured animals.}
    \item{RecapturesSeroplus}{Whether each recaptured animal is seropositive.}
    \item{RecapturesTimeSinceCapt}{For each recaptured animal, the time
      spent since the first capture. }
    \item{RecapturesDurationContamin}{A matrix describing for each
      recaptured animal (rows), whether the year in column belonged to
      the period between the capture and the recapture.}
    \item{BacterialNanalyzed}{The number of analyzed organs for each
      animal in the bacterial cultures dataset.}
    \item{BacterialNpositive}{The number of organs detected as infected
      for each animal in the bacterial cultures dataset.}
    \item{BacterialLogFC}{The logarithm of the CFT titer for each animal
      in the bacterial cultures dataset.}
    \item{BacterialNanimals}{The number of animals in the bacterial
      cultures dataset.}
    \item{BacterialAge}{The age of each animal in the bacterial dataset.}
    \item{CenSeroplusAge}{The age of each animal for which the CFT titer
      was available in the capture dataset.}
    \item{CenSeroplusLogFC}{The log of the CFT titer for each animal for
      which the CFT titer was available in the capture dataset.}
    \item{CenSeroplusYear}{The year for each animal for
      which the CFT titer was available in the capture dataset.}
    \item{CenSeroplusNanimals}{The number of animalsvd for
      which the CFT titer was available in the capture dataset}
    \item{PopSize}{Point estimates of the population size during each
      year of the study.}
    \item{Nyears}{The number of year of the monitoring (6 years in this study).}
  }

}
\details{
  We did not get the authorization to include the original dataset in
  this package, though discussions are still underway to obtain this
  authorization. We therefore simulated this dataset as a replacement,
  hoping that the authorities change their mind soon regarding this
  dataset. 
}
\references{

  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.


}
\examples{

data(BruData)
str(BruData)

}
\keyword{datasets}
