#' print summary for powerNLSEM objects
#' @import stats
#' @importFrom crayon blue
#' @importFrom crayon green
#' @importFrom crayon red
#' @param x Prints the summary of a \code{powerNLSEM} object. \code{x} must be of class \code{"summary.powerNLSEM"}.
#' @param ... Further arguments to use in \code{print}.
#' @return Prints output of summmary of \code{powerNLSEM} object into the console (objects of class \code{summary.powerNLSEM}), but does not change object itself.
#' @exportS3Method print summary.powerNLSEM
#' @export

print.summary.powerNLSEM <- function(x, ...){

     Rok <- ceiling(x$args$R * (1 - x$arg$power_aim)) >= 400

     cat("-----------------------------------------------------------------------------\n")
     cat(paste0(crayon::blue("Model-Implied Simulation-Based Power Estimation: powerNLSEM "),
                packageVersion("powerNLSEM"), "\n\n"))
     cat(crayon::underline("Parameters of Interest (POI):\n"))
     cat(paste(paste0(c(names(x$est)), collapse = ", ")))
     cat(crayon::underline("\n\nTrue Values for POI:\n"))
     print(x$truth)
     cat(crayon::underline("\n\nMethod:\n"))
     cat(paste(x$args$method, collapse = ""))
     cat(crayon::underline("\n\nTest:\n"))
     cat(paste(x$args$test, "z-Test", collapse = ""))
     cat(crayon::underline("\n\nPower (optimized for):\n"))
     cat(paste(ifelse(Rok, crayon::green(x$args$power_aim),
                             crayon::red(x$args$power_aim)), collapse = ""))
     cat(crayon::underline("\n\nType I error/Alpha (for significance decision/z-Tests):\n"))
     cat(paste(x$args$alpha, collapse = ""))
     cat(crayon::underline("\n\nPower Modeling:\n"))
     cat(paste(x$args$power_modeling_method, collapse = ""))
     cat(crayon::underline("\n\nType I error/Alpha (for power modeling):\n"))
     cat(paste(x$args$alpha, collapse = ""))
     cat(crayon::underline("\n\nR (number of replications):\n"))
     cat(paste(ifelse(Rok, crayon::green(x$args$R),
                             crayon::red(x$args$R)), collapse = ""))
     cat(crayon::underline("\n\nConvergence Rate:\n"))
     cat(paste(x$convergenceRate, " (converged samples: ",
               x$args$R*x$convergenceRate,")", collapse = "", sep = ""))
     cat(crayon::underline("\n\nSeed:\n"))
     cat(paste(x$seed, collapse = "", sep = ""))
     cat("\n\n-------------------------------Results---------------------------------------\n")
     cat(crayon::underline("Desired Sample Size:\n"))
     cat(paste(x$N, collapse = ""))
     cat(crayon::underline("\n\nEstimation Performance:\n"))
     print(x$Performance)
     cat(" *weighted bias, absolute bias, relative bias,\n    and root weighted mean squared error\n")
     cat("-----------------------------------------------------------------------------\n")

     return(invisible(x))
}

