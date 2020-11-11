totalPrevalence <-
function(x, limits=c(0.1,0.9), years=2013:2018)
{
    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")

    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame
    nars <- names(rsd)
    if (!any(str_detect(nars,"pit\\["))) {
        stop("Non convenient object")
    }

    mo <- colMeans(rs)
    info <- apply(rs, 2, quantile, limits[1])
    supo <- apply(rs, 2, quantile, limits[2])
    qo <- str_detect(colnames(rs),"pit\\[")

    dgr <- data.frame(Year=years,
                      Prevalence=mo[qo],
                      LowerBound=info[qo],
                      UpperBound=supo[qo])
    return(dgr)
}
