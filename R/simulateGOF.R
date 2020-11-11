simulateGOF <-
function(x, dataset, nsim=sum(sapply(x,nrow)), proba=0.8)
{
    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")

    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame
    na <- names(rsd)

    veri <- c("gammat", "pitu", "lambda", "wt", "sigmaf", "beta0","beta1","beta2","beta3")
    if (!all(sapply(veri, function(x) any(str_detect(na,x)))))
        stop("non-convenient MCMC object")

    verib <- c("CensusesYear", "CensusesNind","CapturesSeroplus", "RecapturesSeroplus",
               "CenSeroplusNanimals", "RecapturesDurationContamin","CenSeroplusAge","CenSeroplusLogFC")
    if (!all(verib%in%names(dataset)))
        stop("non-convenient dataset")

    cll <- (1-proba)/2
    clu <- 1-cll


    ### Check of the number of marked animals
    prm <- rsd %>% select(starts_with("gammat")) %>% as.matrix
    si <- function(i)
    {
        pr <- sapply(1:length(dataset$CensusesYear), function(j) prm[i,dataset$CensusesYear[j]])
        rbinom(length(dataset$CensusesNind), dataset$CensusesNind,pr)
    }
    ss <- sapply(1:nsim, si)
    m <- apply(ss,1,quantile,cll)
    M <- apply(ss,1,quantile,clu)
    Pmarkic <- round(100*mean(dataset$CensusesNmark>=m&dataset$CensusesNmark<=M))


    ## Check of the prevalence of unmarked animals
    pre <- rsd %>% select(starts_with("pitu[")) %>% as.matrix
    si2 <- function(i)
    {
        prec <- sapply(1:length(dataset$CapturesYear), function(j) pre[i,dataset$CapturesYear[j]])
        100*tapply(rbernoulli(length(dataset$CapturesSeroplus), prec), dataset$CapturesYear,mean)
    }
    simm <- sapply(1:nsim, si2)
    obssero <- round(100*tapply(dataset$CapturesSeroplus, dataset$CapturesYear, mean))
    simsero <- apply(simm,1, function(x)
        paste0("[",round(quantile(x,cll)), "-", round(quantile(x,clu)), "]"))


    ## Checke of the prevalence of marked animals
    lambda <- rsd %>% select(starts_with("lambda")) %>% as.matrix
    si3 <- function(i)
    {
        pro <- as.vector(1-exp(-dataset$RecapturesDurationContamin%*%lambda[i,]))
        mean(rbernoulli(length(dataset$RecapturesSeroplus), pro))
    }
    simm3 <- sapply(1:nsim, si3)
    obsserore <- round(100*mean(dataset$RecapturesSeroplus))
    simserore <- paste0("[", round(100*quantile(simm3,cll)),
                        "-", round(100*quantile(simm3,clu)), "]")

    ## Check of the prediction of logFC for the sample of animals where sero was measured
    pba <- rsd %>% select(starts_with("wt")) %>% as.matrix
    si4 <- function(i)
    {
        et <- rbernoulli(dataset$CenSeroplusNanimals, pba[i,][dataset$CenSeroplusYear]) %>% as.numeric
        muelfc <- rsd$beta0[i] + rsd$beta1[i]*et + rsd$beta2[i]*dataset$CenSeroplusAge +
            rsd$beta3[i]*et*dataset$CenSeroplusAge
        rnorm(length(muelfc), mean=muelfc, 1/sqrt(rsd$sigmaf[i]))
    }
    simm3 <- sapply(1:nsim, si4)
    simlogfc <- apply(simm3, 1, quantile, c(cll, clu))
    Plfc <- round(100*mean(dataset$CenSeroplusLogFC>=simlogfc[1,]&dataset$CenSeroplusLogFC<=simlogfc[2,]))

    resu <- list(checkNmarked=Pmarkic, checkSeroCapture=data.frame(obs=obssero, ICsim=simsero),
                 checkSeroRecapture=data.frame(obs=obsserore, ICsim=simserore),
                 checkSeroTiter=Plfc)
    class(resu) <- "gofBIB"
    attr(resu, "nsim") <- nsim
    attr(resu,"proba") <- proba
    return(resu)
}
