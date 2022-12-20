#' Search function to find N for desired power
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Can be LMS or UPI.
#' @param lavModel lavModel object describing the model.
#' @param search_method String stating the search method. "smart" or "bruteforce".
#' @param Ntotal Total number of models to be fitted. Higher number results in higher precision and longer runtime.
#' @param power_aim Minimal power value to approximate. Default to .8.
#' @param alpha Type I-error rate. Default to .05.
#' @param ... Additional arguments passed on to the search functions.
#' @import MplusAutomation
#' @export
#'
power_search <- function(POI,
                         method, lavModel,
                         search_method,
                         Ntotal = 1000,
                         power_aim = .8,
                         alpha = .05, ...)
{
     dotdotdot <- list(...)

     if(!is.null(dotdotdot$N_start)){
          N_start <- dotdotdot$N_start
     }else{
          N_start <- nrow(lavModel)*10
     }
     if(!is.null(dotdotdot$lb)){
          lb <- max(dotdotdot$lb, nro(lavModel))
     }else{
          lb <- nrow(lavModel)
     }
     if(!is.null(dotdotdot$steps)){
          steps <- dotdotdot$steps
     }else{
          steps <- 20
     }
     if(!is.null(dotdotdot$type)){
          type <- dotdotdot$type
     }else{
          type <- "u"
     }

     if(tolower(search_method %in% c("smart", "smart_search")))
     {
          args <- names(formals(smart_search))
          out <- do.call("smart_search", args = mget(args))
     }else if(tolower(search_method %in% c("bruteforce", "brute", "force", "brute force", "brute_force")))
     {
          args <- names(formals(bruteforce_search))
          out <- do.call("bruteforce_search", args = mget(args))
     }

     return(out)
}

