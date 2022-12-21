#' Search function to find N for desired power
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Can be LMS or UPI.
#' @param lavModel lavModel object describing the model.
#' @param search_method String stating the search method. "smart" or "bruteforce".
#' @param Ntotal Total number of models to be fitted. Higher number results in higher precision and longer runtime.
#' @param power_aim Minimal power value to approximate. Default to .8.
#' @param alpha Type I-error rate. Default to .05.
#' @param CORES Number of cores used for parallelization. Default to number of available cores - 2.
#' @param verbose Logical whether progress should be printed in console. Default to TRUE.
#' @param ... Additional arguments passed on to the search functions.
#' @import MplusAutomation
#' @export
#'
power_search <- function(POI,
                         method, lavModel,
                         lavModel_Analysis,
                         data_transformations,
                         search_method,
                         Ntotal = 1000,
                         power_aim = .8,
                         alpha = .05,
                         CORES,
                         verbose,
                         Ns = NULL,
                         N_start = nrow(lavModel)*10,
                         type = "u",
                         steps = 10,
                         lb = nrow(lavModel),
                         switchStep = round(steps/2))
{
     if(tolower(search_method %in% c("smart", "smart_search")))
     {
          args <- names(formals(smart_search))
          args <- args[args != "..."]
          out <- do.call("smart_search", args = mget(args))
     }else if(tolower(search_method %in% c("bruteforce", "brute", "force", "brute force", "brute_force")))
     {
          args <- names(formals(bruteforce_search))
          args <- args[args != "..."]
          out <- do.call("bruteforce_search", args = mget(args))
     }

     return(out)
}

