#' print powerNLSEM objects
#' @param x object of class powerNLSEM
#' @param ... Additional parameters for print
#' @returns \code{powerNLSEM} object
#' @export
#' @exportS3Method


print.powerNLSEM <- function (x, ...)
{
     obj <- unclass(x)
     print(obj[c("N", "alpha", "beta", "power",
                 "convergenceRate",
                 "AveragePerformance", "Performance", "runtime", "call")], ...)
     invisible(x)
}
