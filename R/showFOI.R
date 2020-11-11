showFOI <-
function(x, year0=2012)
{
    if (!inherits(x,"mcmc.list"))
        stop("x should be of class mcmc.list")

    rs <- do.call(rbind, x)
    rsd <- rs %>% as.data.frame()
    na <- names(rsd)

    veri <- c("lambdat")
    if (!all(sapply(veri, function(x) any(str_detect(na,x)))))
        stop("non-convenient MCMC object")

    distriproinf <- rsd %>% as.data.frame() %>% dplyr::select(dplyr::starts_with("lambda")) %>% as.matrix()
    distriproinf <- 1-exp(-distriproinf) %>% as.data.frame()
    names(distriproinf) <- gsub("lambda", "Prop", names(distriproinf))
    toto <- distriproinf %>% tidyr::gather() %>% dplyr::rename(Annee=.data$key,Proportion=.data$value) %>%
        dplyr::mutate(Annee=factor(as.numeric(as.factor(.data$Annee))+year0))

    ggplot2::ggplot(toto,ggplot2::aes(x=factor(.data$Annee), y=.data$Proportion))+
        ggplot2::geom_boxplot(fill="grey")+
        ggplot2::ylab("Prop. newly infected animals")+
        ggplot2::xlab("Year")
}
