#' print powerNLSEM objects
#' @param x object of class powerNLSEM
#' @returns \code{powerNLSEM} object
#' @export
#' @exportS3Method


print.powerNLSEM <- function (x, ...)
{
     obj <- unclass(x)
     print(obj[c(1,3,4,5,6,7,8,9)])
     invisible(x)
}
