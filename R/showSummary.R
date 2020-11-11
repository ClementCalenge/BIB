showSummary <-
function(x, year0=2012)
{

    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")

    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame
    nars <- names(rsd)
    mi <- min(c(length(grep("pit\\[", nars)),
                length(grep("gammat", nars)),
                length(grep("wt", nars))))
    if (mi==0) {
        stop("Non convenient object")
    }

    o <- rsd %>% dplyr::select(starts_with("pit[")) %>% tidyr::gather()
    o$key <- substr(o$key,5,5) %>% as.numeric
    o$key  <- o$key+year0
    p1 <- ggplot2::ggplot(o, ggplot2::aes(x=factor(.data$key), y=.data$value))+
        ggplot2::geom_violin()+
        ggplot2::geom_boxplot(width=0.1, fill="grey")+ggplot2::xlab("Year")+
        ggplot2::ylab("Overall seroprevalence")+
        ggplot2::ggtitle("(A)",
                "")+
        ggplot2::stat_summary(fun=median, colour="red", geom="line", ggplot2::aes(group = 1))


    o <- rsd %>% dplyr::select(starts_with("gammat")) %>% tidyr::gather()
    o$key <- substr(o$key,8,8) %>% as.numeric
    o$key  <- o$key+year0
    p2 <- ggplot2::ggplot(o, ggplot2::aes(x=factor(.data$key), y=.data$value))+ggplot2::geom_violin()+
        ggplot2::geom_boxplot(width=0.1, fill="grey")+ggplot2::xlab("Year")+
        ggplot2::ylab("Prop. marked animals")+
        ggplot2::ggtitle("(B)",
                "")+
        ggplot2::stat_summary(fun=median, colour="red", geom="line", ggplot2::aes(group = 1))


    o <- rsd %>% dplyr::select(starts_with("wt")) %>% tidyr::gather()
    o$key <- substr(o$key,4,4) %>% as.numeric
    o$key  <- o$key+year0
    p3 <- ggplot2::ggplot(o, ggplot2::aes(x=factor(.data$key), y=.data$value))+ggplot2::geom_violin()+
        ggplot2::geom_boxplot(width=0.1, fill="grey")+ggplot2::xlab("Year")+
        ggplot2::ylab("Prop. actively infected sero+")+
        ggplot2::ggtitle("(C)",
                "")+
        ggplot2::stat_summary(fun=median, colour="red", geom="line",
                              ggplot2::aes(group = 1))


    oa <- rsd %>% dplyr::select(starts_with("pit[")) %>% tidyr::gather()
    ob <- rsd %>% dplyr::select(starts_with("wt")) %>% tidyr::gather()
    o <- oa
    o$value <- oa$value*ob$value
    o$key <- substr(o$key,5,5) %>% as.numeric
    o$key  <- o$key+year0
    p4 <- ggplot2::ggplot(o, ggplot2::aes(x=factor(.data$key), y=.data$value))+ggplot2::geom_violin()+
        ggplot2::geom_boxplot(width=0.1, fill="grey")+ggplot2::xlab("Year")+
        ggplot2::ylab("Prop. actively infected")+
        ggplot2::ggtitle("(D)",
                "")+
        ggplot2::stat_summary(fun=median, colour="red", geom="line", ggplot2::aes(group = 1))+
        ggplot2::ylim(c(0,1))


    gridExtra::grid.arrange(p1, p2, p3, p4, nrow=2)
}
