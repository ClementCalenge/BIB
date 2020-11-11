testSpatial <-
function(x, dataSpatial, proba=0.95, year0=2012, nsim=sum(sapply(x,nrow)))
{

    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")
    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame
    na <- names(rsd)
    veri <- "gammat"
    if (!any(str_detect(na,veri)))
        stop("non-convenient MCMC object")

    if (!inherits(dataSpatial, "data.frame"))
        stop("dataSpatial should be a data.frame")

    if (!all(names(dataSpatial)%in%c("Year", "Ntot", "Nmarked", "Zone")))
        stop("dataSpatial does not contain the correct variables")


    cll <- (1-proba)/2
    clu <- 1-cll

    prm <- rsd %>% select(starts_with("gammat")) %>% as.matrix

    sis <- function(i) {
            pr <- sapply(1:nrow(dataSpatial), function(j) prm[i,dataSpatial[j,1]-year0])
            rbinom(nrow(dataSpatial), dataSpatial[,2],pr)
    }

    sss <- sapply(1:nsim, sis)
    sb <- apply(sss,2,function(x) round(100*x/dataSpatial[,2]))
    ms <- apply(sb,1,quantile,cll)
    Ms <- apply(sb,1,quantile,clu)
    ds2 <- as.data.frame(dataSpatial) %>% mutate(ICmin=ms, ICmax=Ms)
    ds2$Nmarked <- paste0(ds2$Nmarked, " (",
                          round(100*dataSpatial[,3]/dataSpatial[,2]),"\\%)")
    ds2$CI <- paste0("[",ds2$ICmin,"\\%--",ds2$ICmax,"\\%]")
    ds2$Zone <- factor(ds2$Zone, labels=c("Southwest","Northeast"))
    ds2$ICmax <- ds2$ICmin <- NULL

    return(ds2)
}
