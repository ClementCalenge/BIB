robustnessPitWit <- function(be1p=1, be2p=1, be1w=1, be2w=1, dataset, niter=5000, thin=5)
{

    mo <- paste0("model {

   ## Priors
   for (i in 1:Nyears) {
     pitu[i]~dbeta(", be1p, ",", be2p,")  ## prevalence among unmarked
     gammat[i]~dbeta(1,1)  ## proportion of marked animals
   }

   ## Force of infection
   eta~dunif(0,100)


   ## Model of infection
   beta0~dnorm(0.0, 0.001)
   beta1~dnorm(0.0, 0.001)
   beta2~dnorm(0.0, 0.001)
   beta3~dnorm(0.0, 0.001)
   sigmaf~dgamma(0.001, 0.001)


   ## Proba to detect infection with bacterioculture
   probad~dbeta(1,1)

   ## Proba of infection by bacteria in the sample
   probai~dbeta(1,1)

   ## Proba of infection per year
   for (i in 1:Nyears) {
      wt[i]~dbeta(",be1w,",", be2w, ")
   }

   ############
   ## Model relating logFC et infection
   for (i in 1:BacterialNanimals) {
      probaDetection[i] <- probad*e[i]
      BacterialNpositive[i]~dbin(probaDetection[i], BacterialNanalyzed[i])
      e[i]~dbern(probai)
      mulfc[i] <- beta0 + beta1*e[i] + beta2*BacterialAge[i] + beta3*e[i]*BacterialAge[i]
      BacterialLogFC[i]~dnorm(mulfc[i], sigmaf)
   }

   ############
   ## Model using the coefficient used for the proportion of infection
   for (i in 1:CenSeroplusNanimals) {
     muelfc[i] <- beta0 + beta1*et[i] + beta2*CenSeroplusAge[i] + beta3*et[i]*CenSeroplusAge[i]
     CenSeroplusLogFC[i]~dnorm(muelfc[i], sigmaf)
     et[i]~dbern(wt[CenSeroplusYear[i]])
   }


   ## For each counting day
   for (i in 1:CensusesNdays) {
      CensusesNmark[i]~dbin(gammat[CensusesYear[i]], CensusesNind[i]) ## number of marked animals
   }

   ## Inference on the prevalence of the captures and recaptures
   for (i in 1:CapturesN) {
      CapturesSeroplus[i]~dbern(pitu[CapturesYear[i]])
   }


   ## Estimation sero recapture
   expos <- RecapturesDurationContamin%*%lambdat
   for (i in 1:RecapturesN) {
      prevage[i] <- 1-exp(-expos[i])
      RecapturesSeroplus[i]~dbern(prevage[i])
   }

   ## Calculation of prevalence
   pitm[1] <- 0
   Nrsurviv[1] <- CaptureNumber[1]  ## Marked of the year
   Nrinfect[1] <- 0
   ## For each additional year
   for (i in 2:Nyears) {
      whichsurv[i-1] ~ dbin(0.95, Nrsurviv[i-1])
      Nrsurviv[i] <- whichsurv[i-1]+CaptureNumber[i] ## The animals that survivedd from last year +
                                                ## marked of the year
                                                ## (because first year, marked
                                                ## available to observation)
      Nrinfect[i] ~ dbin(1-exp(-lambdat[i]), whichsurv[i-1])
      pitm[i] <- Nrinfect[i]/Nrsurviv[i]
   }

   ## Inference on total prevalence
   for (i in 1:Nyears) {
       pit[i] <- ((1-gammat[i]) * pitu[i])+(gammat[i] * pitm[i])
   }

   ## Force of infection of marked animals
   for (i in 1:Nyears) {
     lambdat[i] <- eta*ut[i]
     ut[i] <- (1-gammat[i])*pitu[i]*wt[i]
   }

   ## Proportion of the population of females with active infection
   for (i in 1:Nyears) {
      omegat[i] <- pit[i]*wt[i]
   }



}"
)

    obmcy <- jags.model(textConnection(mo),
                        data=dataset, n.chains = 5,
                        inits =
                            list(e=as.numeric(dataset$BacterialNpositive>0),
                                 eta=0.1))

    sam <- coda.samples(obmcy, n.iter=niter, thin=thin,
                        variable.names = c("pitu","pitm","pit", "gammat",
                                           "eta","beta0","beta1","beta2","beta3",
                                           "sigmaf","probad","probai",
                                           "wt","lambdat", "omegat"))

    return(sam)
}
