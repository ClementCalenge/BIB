getTable2and3 <-
function(x, years=2013:2018, proba=0.9, table2=TRUE)
{
    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")


    cll <- (1-proba)/2
    clu <- 1-cll

    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame
    na <- names(rsd)

    veri <- c("pitu", "pitm", "pit\\[", "gammat")
    if (!all(sapply(veri, function(x) any(str_detect(na,x)))))
        stop("non-convenient MCMC object")


    foob <- function(nab)
    {
        na <- paste0(nab,"\\[")
        foo <- function(x) c(mean(x), quantile(x, c(cll,clu)))
        aa <- apply(rsd,2,foo) %>% t %>% as.data.frame %>% rownames_to_column
        uu <- aa %>% filter(str_detect(aa$rowname, na))
        names(uu) <- c("name","val","lb","ub")
        uu <- uu %>% dplyr::mutate(Year=years) %>%
            dplyr::mutate(res=str_c(round(.data$val,2), " [",
                                    round(.data$lb,2), "-",
                                    round(.data$ub,2), "]")) %>%
            dplyr::select(.data$Year, .data$res)
        return(uu)
    }


    if (table2) {
        prt <- foob("pit")
        prm <- foob("pitu")
        prr <- foob("pitm")
        prom <- foob("gammat")
        ddd <- data.frame(Year=prt$Year, Overall.Prevalence=prt$res,
                          Seroprev.marked=prr$res, Seroprev.unmarked=prm$res,
                          Prop.Marked=prom$res)
        ddd1 <- ddd
    } else {
        prb <- foob("wt")

        oa <- rsd %>% dplyr::select(starts_with("pit[")) %>% as.matrix
        ob <- rsd %>% dplyr::select(starts_with("wt")) %>% as.matrix
        ppor <- oa*ob %>% as.data.frame
        names(ppor) <- paste0("ppor[",1:ncol(ppor),"]")
        rsd <- cbind(rsd, ppor)
        pro <- foob("ppor")
        ddd <- data.frame(Year=prb$Year, Prop.sero.active.inf=prb$res,
                          Overall.prop.active.inf=pro$res)
    }
    return(ddd)
}
