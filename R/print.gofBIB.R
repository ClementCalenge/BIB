print.gofBIB <-
function(x, ...)
{
    if (!inherits(x, "gofBIB"))
        stop("x should inherit the class \"gofBIB\"")

    cat("######################################\n##\n## Goodness of fit of the BIB model\n\n")
    cat("Number of simulations:", attr(x, "nsim"),"\n")
    cat("Probability used for the credible intervals in this check:", 100*attr(x, "proba"),"%\n\n")
    cat("Percentage of the simulate CI including\nthe obs. value of number of marked animals:", x$checkNmarked,"%\n\n")

    cat("Comparison of the observed percentage of seropositive animals among\ncaptured unmarked animals with the simulated CI for each year:\n")
    print(x$checkSeroCapture)
    cat("\n\n")


    cat("Comparison of the observed percentage of seropositive animals among\ncaptured marked animals with the simulated CI (all years pooled):\n")
    print(x$checkSeroRecapture)
    cat("\n\n")

    cat("Percentage of logFC values in the simulated CI for analysed animals:", x$checkSeroTiter, "%\n\n")

}
