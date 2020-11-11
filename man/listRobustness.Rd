\name{listRobustness}
\alias{listRobustness}
\docType{data}
\title{
  Test of the Robustness of our Model to the Changes in Prior Distribution
}
\description{
  This list contains the result of the simulations carried out to test
  the sensitivity of our fit to changes in the prior for the
  seroprevalence of unmarked animals and for the prevalence of active
  infection among seropositive animals (see example for the code used to
  create this list).
}
\usage{data("listRobustness")}
\format{

  This is a list of three data.frames containing: (i) the estimated
  total seroprevalence (and lower/upper bounds of the 80\% credible
  interval on this prevalence) of Brucellosis when the shape parameters
  of the beta prior distribution on the seroprevalence of unmarked
  animals varies from 1 to 5; (ii) the estimated proportion of actively
  infected animals among seropositive animals (and lower/upper bounds of
  the 80\% credible interval on this prevalence) of Brucellosis when the
  shape parameters of the beta prior distribution on this proportion
  varies from 1 to 5; (iii) the estimated proportion of actively
  infected animals among all animals (and lower/upper bounds of the 80\%
  credible interval on this prevalence) of Brucellosis when the shape
  parameters of the two beta prior distributions on both the proportion
  of seropositive animals that are actively infected, and the
  seroprevalence of unmarked animals varies from 1 to 5.
  
}
\references{

  Calenge C., Lambert S., Petit E., Thebault A., Gilot-Fromont E., Toigo
  C., Rossi S. (in prep.) Estimating disease prevalence and temporal
  dynamics using biased capture serological data in a wildlife
  reservoir: the example of brucellosis in Alpine ibex (Capra
  ibex). Submitted to Preventive Veterinary Medicine.

}
\examples{

\dontrun{

## The set of shape parameters
listebeta <- c(1,1.5,2,5)

#######
## Total seroprevalence = f(prior seroprev. among unmarked)

so <- lapply(listebeta, function(x) {
    samb <- robustnessPitWit(x,x,1,1, BruData)
    totalPrevalence(samb)
})
dfso <- do.call(rbind, so)
dfso$beta <- rep(listebeta,each=6)


#######
## Prop. Active infect. among sero+ = f(prior Prop. Active infect. among sero+)

sow <- lapply(listebeta, function(x) {
    samb <- robustnessPitWit(1,1,x,x, BruData)
    rs <- do.call(rbind,samb)
    r1 <- rs[,str_detect(colnames(rs), "wt")]
    data.frame(Year=2013:2018, wtm=colMeans(r1),
               lb=apply(r1,2,quantile,c(0.1)),
               ub=apply(r1,2,quantile,c(0.9)))
})
dfsow <- do.call(rbind, sow)
dfsow$beta <- rep(listebeta,each=6)


#######
## Prop. Active infect. in pop = f(prior Prop. Active infect. among
## sero+, prior seroprev among unmarked)

soom <- lapply(listebeta, function(x) {
    samb <- robustnessPitWit(x,x,x,x, BruData)
    rs <- do.call(rbind,samb)
    r1a <- rs[,str_detect(colnames(rs), "wt")]
    r1b <- rs[,str_detect(colnames(rs), "pit\\[")]
    r1 <- r1a*r1b
    data.frame(Year=2013:2018, omegam=colMeans(r1),
               lb=apply(r1,2,quantile,c(0.1)), ub=apply(r1,2,quantile,c(0.9)))
})
dfsoom <- do.call(rbind, soom)
dfsoom$beta <- rep(listebeta,each=6)



listRobustness <- list(pit=dfso, wit=dfsow, pitwit=dfsoom)
}

}
\keyword{datasets}
