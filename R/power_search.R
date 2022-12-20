#' Search function to find N for desired power
#' @param method Method used to fit to the data.
#' @param lavModel lavModel object describing the model.
#' @import MplusAutomation
#' @export
power_search <- function(method, lavModel,
                         search_method,
                         N_start = 200, type = "u",
                         Ntotal = 1000, steps = 10,
                         power_aim = .8, alpha = .05, ...)
{
     dotdotdot <- list(...)
}



test <- function(x, ...)
{
     dotdotdot <- list(...)
     print(dotdotdot)
     return(x+1)
}
test(x = 2, a = 21)
