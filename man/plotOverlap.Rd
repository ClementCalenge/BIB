\name{plotOverlap}
\alias{plotOverlap}
\title{
  Table Summarizing the Overlap Between Prior and Posterior Distributions
}
\description{
  This function reproduces the table presented in the appendix of the
  paper of Calenge et al. (in prep) showing the table summarizing the
  overlap between the prior and the posterior distributions of the
  parameters. 
}
\usage{
plotOverlap(x, df, cex = 1.5)
}
\arguments{
  \item{x}{
    An object of class \code{mcmc.list} containing the same parameters
    as in \code{df}.
  }
  \item{df}{
    A data.frame (see details).
  }
  \item{cex}{
    the size of the text in the table
  }
}
\details{
  
  The object \code{df} must contain the following variables: (i)
  \code{Parameter}: the name of the parameter in the model \code{i},
  (ii) \code{tau}: the value of the overlap parameter tau between the
  prior and the posterior,(iii) \code{namo}: a character string
  containing a plotmath expression to be used in the first column of the
  resulting table, (iv) \code{prior}: a character string indicating the
  form of the prior for the parameter; note that the parameter should
  be named \code{z} in this expression, (v) three columns named
  \code{from,to,n} passed to the function \code{density}, and used to
  calculate the values of the parameter over which the prior and
  posterior should be compared.
  
}
\value{
  Nothing.  Just a plot.
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
  \code{\link{overlapPriorPost}} to calculate the coefficient
  \code{tau}.
}
\examples{

## Use the MCMC samples from the model ModelCatalyticYear
rs <- do.call(rbind, SamplesModelCatalyticYear)

## Detects the name of the parameters for which tau should be calculated
library(stringr)
nam <- c("pitu","gammat","probad","wt","beta","^eta","sigmaf")
nam <- unlist(lapply(nam, function(z) colnames(rs)[str_detect(colnames(rs),z)]))


## For each parameter, calculates the value of tau:

lb <- sapply(nam, function(na) {
    if (any(sapply(c("pitu","gammat","probad","probai","wt"), function(x) str_detect(na, x))))
        tau <- overlapPriorPost(SamplesModelCatalyticYear, na, prior=1, from=0, to=1)
    if (str_detect(na, "beta"))
        tau <- overlapPriorPost(SamplesModelCatalyticYear, na)
    if (str_detect(na, "eta"))
        tau <- overlapPriorPost(SamplesModelCatalyticYear, na, prior=1/100, from=0, to=100)
    if (str_detect(na, "sigmaf"))
        tau <- overlapPriorPost(SamplesModelCatalyticYear, na,
                                prior="dgamma(den$x[j],0.001, 0.001)", from=0, to=1000)
    return(tau)
})


## Builds the data.frame required by plotOverlap
df <- data.frame(Parameter=nam, tau=lb)
namo <- paste0("expression(",c(paste0("pi[",1:6,"]^(u)"), paste0("gamma[",1:6,"]"), "p[d]",
                               paste0("w[",1:6,"]"), paste0("beta[",0:3,"]"),
                               "eta","sigma[f]"),")")
df$namo <- namo
df$prior <- "dunif(z)"
df$prior[20:23] <- "dnorm(z,0,sqrt(1000))"
df$prior[24] <- "dunif(z,0,100)"
df$prior[25] <- "dgamma(z,0.001,scale=0.001)"
df$from <- 0
df$to <- 1
df$from[20:23] <- -4
df$to[20:23] <- 4
df$from[24] <- 0
df$to[24] <- 30
df$from[25] <- 0.01
df$to[25] <- 5


plotOverlap(SamplesModelCatalyticYear, df, cex=1)



}
\keyword{ hplot }
