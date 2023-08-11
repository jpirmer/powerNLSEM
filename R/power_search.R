#' Search function to find N for desired power
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Can be LMS or UPI.
#' @param lavModel lavModel object describing the model.
#' @param lavModel_Analysis lavModel object containg the parameters to be estimated.
#' @param data_transformations Object containing info on data transformations.
#' @param power_modeling_method Power modeling method used to model significant parameter estimates. Default to \code{"probit"} indicating glm with probit link function with sqrt(n) as predictor. Alternative is \code{"logit"}.
#' @param search_method String stating the search method. Default to \code{"smart"}. Alternative is \code{"bruteforce"}.
#' @param R Total number of models to be fitted. Higher number results in higher precision and longer runtime.
#' @param power_aim Minimal power value to approximate. Default to .8.
#' @param alpha Type I-error rate. Default to .05.
#' @param CORES Number of cores used for parallelization. Default to number of available cores - 2.
#' @param verbose Logical whether progress should be printed in console. Default to TRUE.
#' @param  Ns Sample sizes used in power estimation process. Default to \code{NULL}.
#' @param N_start Starting sample size. Default to \code{nrow(lavModel)*10}
#' @param type  Default to \code{"u"}.
#' @param steps Steps used in search_method = "smart", i.e., the smart algorithm. This is ignored if bruteforce is used. Default to 10.
#' @param lb Lower bound of N used in search. Default to \code{nrow(lavModel)}
#' @param switchStep Steps after which smart search method changes from exploration to exploitation. Default to \code{round(steps/2)}
#' @param uncertainty_method Uncertainty method used for confidence intervals. Default to \code{""}
#' @param seeds Seeds for reproducibility.
#' @export

power_search <- function(POI,
                         method, lavModel,
                         lavModel_Analysis,
                         data_transformations,
                         search_method,
                         power_modeling_method,
                         R = 1000,
                         power_aim = .8,
                         alpha = .05,
                         CORES,
                         verbose,
                         Ns = NULL,
                         N_start = nrow(lavModel)*10,
                         type = "u",
                         steps = 10,
                         lb = nrow(lavModel),
                         switchStep = round(steps/2),
                         uncertainty_method = "",
                         seeds)
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

