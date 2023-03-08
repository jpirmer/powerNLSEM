#' powerNLSEM function
#' @param model Model in lavaan syntax. See documentation for help and examples.
#' @param POI Parameter Of Interest as a vector of strings. Must be in lavaan-syntax without any spaces. Nonlinear effects should have the same ordering as in model.
#' @param method Method used to fit to the data. Can be LMS or UPI.
#' @param power_modeling_method Power modeling method used to model significant parameter estimates. Default to \code{"logit"} indicating glm with logit link function with sqrt(n).
#' @param search_method String stating the search method. Can be \code{"smart"} or \code{"bruteforce"}.
#' @param Ntotal Total number of models to be fitted. Higher number results in higher precision and longer runtime. Default to 2000.
#' @param power_aim Minimal power value to approximate. Default to .8.
#' @param alpha Type I-error rate. Default to .05.
#' @param CORES Number of cores used for parallelization. Default to number of available cores - 2.
#' @param verbose Logical whether progress should be printed in console. Default to TRUE.
#' @param seed Seed for replicability. Default to \code{NULL}, then a seed is drawn at random, which will also be saved in the output.
#' @param ... Additional arguments passed on to the search functions.
#' @export

powerNLSEM <- function(model, POI,
                       method,
                       power_modeling_method = "logit",
                       search_method,
                       Ntotal = 2000,
                       power_aim = .8,
                       alpha = .05,
                       CORES = max(c(parallel::detectCores()-2, 1)),
                       verbose = T, seed = NULL,
                       ...)
{
     if(is.null(seed)) seed <- sample(1:10^9, size = 1)
     set.seed(seed)
     t0 <- proc.time()

     ### prepare model ----
     lavModel <- lavaan:::lavMatrixRepresentation(lavaan::lavaanify(model = model,
                                                                    meanstructure = T,auto.var = T,
                                                                    auto.cov.lv.x = T, auto.cov.y = T,
                                                                    as.data.frame. = T))
     lavModel <- add_varType(lavModel)
     temp <- check_model(lavModel); lavModel <- temp$lavModel; added_model_syntax <- temp$added_model_syntax
     Manifests <- handle_manifests(lavModel = lavModel, treat_manifest_as_latent = "ov")
     lavModel_Analysis <- Manifests$lavModel_Analysis
     data_transformations <- Manifests$data_transformations
     lavModel_Analysis$matchLabel <- toupper(paste0(lavModel_Analysis$lhs,
                                                    lavModel_Analysis$op,
                                                    lavModel_Analysis$rhs))

     ### initialize ----
     dotdotdot <- list(...)
     if(!is.null(dotdotdot$Ns)){
          Ns <- dotdotdot$Ns
          search_method <- "bruteforce"
     }else{
          Ns <- NULL
     }
     dotdotdot <- list(...)

     if(!is.null(dotdotdot$N_start)){
          N_start <- dotdotdot$N_start
     }else{
          N_start <- nrow(lavModel)*10
     }
     if(!is.null(dotdotdot$lb)){
          lb <- max(dotdotdot$lb, nrow(lavModel))
     }else{
          lb <- nrow(lavModel)
     }
     if(!is.null(dotdotdot$steps)){
          steps <- dotdotdot$steps
     }else{
          steps <- 10
     }
     if(!is.null(dotdotdot$switchStep)){
          switchStep <- dotdotdot$switchStep
     }else{
          switchStep <- round(steps/2)
     }
     if(!is.null(dotdotdot$type)){
          type <- dotdotdot$type
     }else{
          type <- "u"
     }
     if(!is.null(dotdotdot$uncertainty_method)){
          uncertainty_method <- dotdotdot$uncertainty_method
     }else{
          uncertainty_method <- ""
     }

     POI <- stringr::str_replace_all(string = POI, pattern = " ", replacement = "")


     args <- names(formals(power_search))
     args <- args[args!="..."]


     if(!all(POI %in% lavModel_Analysis$matchLabel)){
          POI_missing <- POI[(!POI %in% lavModel_Analysis$matchLabel)]
          stop(paste0("Parameters Of Interest (POI) not given in model. Includes: ",
                      paste(POI_missing, collapse = ", ", sep = ""), "!"))
     }

     ### run power analysis ----
     out <- do.call("power_search", mget(args))


     ### return results ----
     t <- proc.time()-t0
     out$power <- power_aim
     out$beta <- 1-power_aim
     out$alpha <- alpha
     out$search_method <- search_method
     out$power_modeling_method <- power_modeling_method
     out$runtime <- t
     out$seed <- seed
     out$model <- paste0(model, added_model_syntax, collapse = "\n")
     class(out) <- c("powerNLSEM", "list")
     return(out)
}
