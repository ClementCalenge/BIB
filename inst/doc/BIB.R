## ----setup, include=FALSE, cache=FALSE--------------------
# set global chunk options
library('knitr')
opts_chunk$set(fig.path="caperpy-",
               fig.align="center",
               fig.show="hold",
               echo=TRUE,
               results="markup",
               fig.width=10,
               fig.height=10, out.width='\\linewidth',
               out.height='\\linewidth',
               cache=FALSE,
               dev='png',
               concordance=TRUE,
               error=FALSE)
opts_knit$set(aliases = c(h = 'fig.height',
              w = 'fig.width',
              wo='out.width',
              ho='out.height'))
options(replace.assign=TRUE,width=60)
set.seed(9567)

## ----eval=FALSE-------------------------------------------
#  ## If devtools is not yet installed, type
#  install.packages("devtools")
#  
#  ## Install the package caperpyogm
#  devtools::install_github("ClementCalenge/BIB")

## ----load-library-----------------------------------------
library(BIB)

## ----description-lekcounts--------------------------------
str(BruData)

## ----JAGS-model-paper, eval=FALSE-------------------------
#  cat(ModelCatalyticYear)

## ----eval=FALSE-------------------------------------------
#  library(rjags)
#  
#  obmcy <- jags.model(textConnection(ModelCatalyticYear),
#                      data=BruData, n.chains = 5,
#                      inits =
#                          list(e=as.numeric(BruData$BacterialNpositive>0),
#                               eta=0.1))
#  
#  SamplesModelCatalyticYear <-
#      coda.samples(obmcy, n.iter=5000, thin=5,
#                   variable.names = c("pitu","pitm","pit", "gammat",
#                                      "eta","beta0","beta1","beta2","beta3",
#                                      "sigmaf","probad","probai",
#                                      "wt","lambdat", "omegat"))
#  

## ----traceplot--------------------------------------------
library(bayesplot)
library(MCMCvis)
cm <- MCMCchains(SamplesModelCatalyticYear,
                 params=c("pitu","pitm","pit", "gammat",
                          "eta","beta0","beta1","beta2","beta3",
                          "sigmaf","probad","probai",
                          "wt","lambdat", "omegat"),
                 mcmc.list=TRUE)
mcmc_trace(cm)

## ----gelman-diagnostic, eval=FALSE------------------------
#  gelman.diag(SamplesModelCatalyticYear, multivariate=FALSE)

## ----simulate-datasets------------------------------------
simg <- simulateGOF(SamplesModelCatalyticYear, BruData)
simg

## ----plot-prev--------------------------------------------
library(ggplot2)
df <- do.call(rbind, lapply(list(SamplesModelCatalyticYear,
                                 SamplesModelConstPrev,
                                 SamplesModelCatalytic,
                                 SamplesModelDensity),
                            totalPrevalence))

df$Model <- factor(c(rep("Variable FOI (used in the paper)", 6),
                                 rep("Constant seroprevalence", 6),
                                 rep("Constant force infection", 6),
                                 rep("Variable FOI (approx. density)", 6)),
                                levels=c("Constant seroprevalence",
                                         "Constant force infection",
                                         "Variable FOI (used in the paper)",
                                         "Variable FOI (approx. density)"))

pd <- position_dodge(0.2)
ggplot(df,aes(x=Year, y=Prevalence, col=Model, group=Model))+
     geom_point(size=3, position=pd)+geom_line(position=pd)+
     geom_errorbar(aes(ymin=LowerBound, ymax=UpperBound), width=0.2, position=pd)+
     xlab("Year")+ylab("Overall seroprevalence")


## ---------------------------------------------------------
getTable2and3(SamplesModelCatalyticYear)

## ---------------------------------------------------------
getTable2and3(SamplesModelCatalyticYear, table2=FALSE)

## ---------------------------------------------------------
showSummary(SamplesModelCatalyticYear)

## ---------------------------------------------------------
showFOI(SamplesModelCatalyticYear)

## ---------------------------------------------------------
capturesSpatial

## ---------------------------------------------------------
testSpatial(SamplesModelCatalyticYear, capturesSpatial)

## ---------------------------------------------------------
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
namo <- paste0("expression(",c(paste0("pi[",1:6,"]^(u)"), paste0("gamma[",1:6,"]"), "p[d]", paste0("w[",1:6,"]"), paste0("beta[",0:3,"]"), "eta","sigma[f]"),")")
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

## ----eval=FALSE-------------------------------------------
#  library(stringr)
#  
#  listebeta <- c(1,1.5,2,5)
#  
#  so <- lapply(listebeta, function(x) {
#      samb <- robustnessPitWit(x,x,1,1, BruData)
#      totalPrevalence(samb)
#  })
#  dfso <- do.call(rbind, so)
#  dfso$beta <- rep(listebeta,each=6)
#  sow <- lapply(listebeta, function(x) {
#      samb <- robustnessPitWit(1,1,x,x, BruData)
#      rs <- do.call(rbind,samb)
#      r1 <- rs[,str_detect(colnames(rs), "wt")]
#      data.frame(Year=2013:2018, wtm=colMeans(r1),
#                 lb=apply(r1,2,quantile,c(0.1)),
#                 ub=apply(r1,2,quantile,c(0.9)))
#  })
#  dfsow <- do.call(rbind, sow)
#  dfsow$beta <- rep(listebeta,each=6)
#  soom <- lapply(listebeta, function(x) {
#      samb <- robustnessPitWit(x,x,x,x, BruData)
#      rs <- do.call(rbind,samb)
#      r1a <- rs[,str_detect(colnames(rs), "wt")]
#      r1b <- rs[,str_detect(colnames(rs), "pit\\[")]
#      r1 <- r1a*r1b
#      data.frame(Year=2013:2018, omegam=colMeans(r1),
#                 lb=apply(r1,2,quantile,c(0.1)), ub=apply(r1,2,quantile,c(0.9)))
#  })
#  dfsoom <- do.call(rbind, soom)
#  dfsoom$beta <- rep(listebeta,each=6)
#  
#  listRobustness <- list(pit=dfso, wit=dfsow, pitwit=dfsoom)

## ---------------------------------------------------------
library(gridExtra)
pd <- position_dodge(0.2)
pl1 <- ggplot(listRobustness$pit,aes(x=Year, y=Prevalence, col=beta, group=beta))+
    geom_point(size=3, position=pd)+geom_line(position=pd)+
    geom_errorbar(aes(ymin=LowerBound, ymax=UpperBound), width=0.2, position=pd)+
    xlab("Year")+ylab("Overall seroprevalence")+ggtitle("(A)")+
    labs(col = "Shape")

pd <- position_dodge(0.2)
pl2 <- ggplot(listRobustness$wit,aes(x=Year, y=wtm, col=beta, group=beta))+
    geom_point(size=3, position=pd)+geom_line(position=pd)+
    geom_errorbar(aes(ymin=lb, ymax=ub), width=0.2, position=pd)+
    xlab("Year")+ylab("Prop. active infection among sero+")+
    ggtitle("(B)")+labs(col = "Shape")

pd <- position_dodge(0.2)
pl3 <- ggplot(listRobustness$pitwit,aes(x=Year, y=omegam, col=beta, group=beta))+
    geom_point(size=3, position=pd)+geom_line(position=pd)+
    geom_errorbar(aes(ymin=lb, ymax=ub), width=0.2, position=pd)+
    xlab("Year")+ylab("Prop. active Infection")+
    ggtitle("(C)")+labs(col = "Shape")

grid.arrange(pl1, pl2, pl3, nrow=2)

