#' Search function to find N for desired power
#' @description
#' The function that initializes the search process. The \code{powerNLSEM} function actually is a wrapper function for \code{power_search}.
#'
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Can be LMS or UPI.
#' @param test Should the parameter be tested with a directed hypothesis (onesided) or with an undirected hypothesis (twosided, also equivalent to Wald-Test for single parameter). Default to \code{"onesided"}.
#' @param lavModel lavModel object describing the model.
#' @param lavModel_Analysis lavModel object containg the parameters to be estimated.
#' @param data_transformations Object containing info on data transformations.
#' @param power_modeling_method Power modeling method used to model significant parameter estimates. Default to \code{"probit"} indicating glm with probit link function with sqrt(n) as predictor. Alternative is \code{"logit"}.
#' @param search_method String stating the search method. Default to \code{"adaptive"} (synonyme is \code{"smart"}). Alternative is \code{"bruteforce"}.
#' @param R Total number of models to be fitted. Higher number results in higher precision and longer runtime.
#' @param power_aim Minimal power value to approximate. Default to .8.
#' @param alpha Type I-error rate for significance decision. Default to \code{.05}.
#' @param alpha_power_modeling Type I-error rate for confidence band around predicted power rate. Used to ensure that the computed \code{N} keeps the desired power value (with the given Type I-error rate \code{alpha_power_modeling} divided by 2). If set to 1, no confidence band is used. Default to \code{.05}.
#' @param CORES Number of cores used for parallelization. Default to number of available cores - 2.
#' @param verbose Logical whether progress should be printed in console. Default to TRUE.
#' @param  Ns Sample sizes used in power estimation process. Default to \code{NULL}.
#' @param N_start Starting sample size for smart algorithm. Default to  \code{10*nrow(lavModel[lavModel$op != "~1", ])} (10 times the number of parameters, excluding the mean structure, without the generation of e.g., factor scores or product indicators).
#' @param distRj  Indicator how the samples sizes should be used in the steps of the smart algorithm: \code{"u"} for many to few to many, \code{"increasing"} for increasing replications and \code{"even"} for evenly distributed replications across steps. Default to \code{"u"}.
#' @param steps Steps used in \code{search_method = "smart"}, i.e., the smart algorithm. This is ignored if bruteforce is used. Default to 10.
#' @param nlb Lower bound of N used in search. Default to \code{5*nrow(lavModel[lavModel$op != "~1", ])} (5 times the number of parameters, excluding the mean structure, in the model without the generation of e.g., factor scores or product indicators), however, some methods can deal with much smaller sample sizes so this can be adjusted. The rule of thumb of 5 times number of parameters is motivated by Wolf et al. (2013)
#' @param switchStep Steps after which smart search method changes from exploration to exploitation. Default to \code{round(steps/2)}. Exploration phase searches for the interval for N so that the resulting power is within \code{[.15, .85]} since the power curve is steepest at .5 and becomes less step towards plus/min \code{Inf}. Exploitation phase searches for an interval for N around the \code{power_aim} argument which shrinks from plus/minus .1 to .01. If \code{swicthStep = Inf}, then only exploration is used. If \code{switchStep} is used then the search process is reset at that point, which results in a new estimation in the bounds of the interval of N independent of the previous ones which might be restricted in change (see also argument( \code{constrainRelChange}).
#' @param constrainRelChange Logical whether the change in the bounds of the interval for N using the smart algorithm should be constrained. This prevents divergence (which is especially an issue for small effect sizes and small \code{R}) but results in biased estimates if the number of steps is too small. Default to \code{TRUE}.
#' @param matchPI Logical passed to \code{semTools::indProd} in order to compute the product indicators: Specify TRUE to use match-paired approach (Marsh, Wen, & Hau, 2004). If FALSE, the resulting products are all possible products. Default to \code{TRUE}. The observations are matched by order given when specifying the measurement model.
#' @param PIcentering String indicating which method of centering should be used when constructing product indicators. String is converted to the arguments \code{meanC}, \code{doubleMC}, and \code{residualMC}, of the \code{semTools::indProd} function. Default to \code{"doubleMC"} for double mean centering the resulting products (Lin et. al., 2010). Use \code{"meanC"} for mean centering the main effect indicator before making the products or \code{"residualC"} for residual centering the products by the main effect indicators (Little, Bovaird, & Widaman, 2006). \code{"none"} or any other input than the previously described results in no centering (use with caution!).
#' @param liberalInspection Logical whether the inspection of estimation truthworthiness should be very liberal (i.e., allowing for non-positive definite Hessians in standard error estimation or non-positive residual covariance matrices or latent covariance matrices). Default to \code{FALSE}. Being liberal is not adviced and should be checked for a single data set!
#' @param FSmethod Method to be used to extract factor scores. Default to \code{"SL"} for the Skrondal and Laake approach that uses regression (\code{"regression"}) factor scores for the independendent variables and \code{"Bartlett"} factor scores for the dependent variables.
#' @param seeds Seeds for reproducibility.
#' @param pathLMS path where (temporal) data and scripts for running LMS using Mplus are stored (using \code{MplusAutomation}). Default to \code{NULL}, then \code{tempdir()} is used.
#' @return Returns a \code{list} that includes the results on model-implied simulation-based power estimation.
#' @references Wolf, E. J., Harrington, K. M., Clark, S. L., & Miller, M. W. (2013). Sample Size Requirements for Structural Equation Models: An Evaluation of Power, Bias, and Solution Propriety. _Educational and Psychological Measurement, 76_(6), 913–934. \doi{10.1177/0013164413495237}
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
                         alpha_power_modeling = .05,
                         CORES,
                         verbose,
                         Ns = NULL,
                         N_start = nrow(lavModel[lavModel$op != "~1", ])*10,
                         distRj = "increasing",
                         steps = 10,
                         nlb = nrow(lavModel[lavModel$op != "~1", ])*5,
                         switchStep = round(steps/2),
                         FSmethod = "SL",
                         test = "onesided",
                         matchPI =TRUE,
                         PIcentering = "doubleMC",
                         liberalInspection = FALSE,
                         constrainRelChange = TRUE,
                         seeds,
                         pathLMS = tempdir())
{
     if(tolower(search_method %in% c("smart", "adaptive")))
     {
          args <- names(formals(smart_search))
          args <- args[args != "..."]
          out <- do.call("smart_search", args = mget(args))
     }else if(tolower(search_method %in% c("bruteforce")))
     {
          args <- names(formals(bruteforce_search))
          args <- args[args != "..."]
          out <- do.call("bruteforce_search", args = mget(args))
     }

     return(out)
}

