#' @import stats

sim_and_fit <- function(n, POI, method,
                        alpha,
                        lavModel, lavModel_Analysis,
                        lavModel_attributes, matrices, data_transformations,
                        prefix,
                        sim_seed,
                        ...){
     set.seed(sim_seed); df_POI <- data.frame("matchLabel" = POI)
     data <- simulateNLSEM(n = n, lavModel = lavModel,
                           lavModel_attributes = lavModel_attributes,
                           matrices = matrices)
     if(tolower(method) == "lms")
     {
          # Mplus only works with upper cases: overwrite df_POI
          df_POI <- data.frame("matchLabel" = toupper(POI))
          lavModel_Analysis$matchLabel <- toupper(lavModel_Analysis$matchLabel)
          fit <- try(LMS(lavModel_Analysis = lavModel_Analysis, data = data,
                         data_transformations = data_transformations, prefix = prefix), silent = TRUE)
     }else if(tolower(method) %in% c("path", "regression", "scaleregression", "sr", "reg")){
          fit <- try(SR(lavModel_Analysis = lavModel_Analysis, data = data,
                        data_transformations = data_transformations), silent = TRUE)
     }
     if(!inherits(fit, "try-error"))
     {
          fit_temp <- fit[fit$matchLabel %in% POI,,drop = FALSE]
          # sort by POI
          fit <- merge(df_POI, fit_temp, by.x = "matchLabel", sort = FALSE)
          out <- try(abs(fit$est/fit$se), silent = TRUE)
          if(!inherits(out, "try-error"))
          {
               out <- out > qnorm(p = 1-alpha)
          }else{
               out <- rep(NA, length(POI))
          }
     }else{
          out <- rep(NA, length(POI))
     }
     names(out) <- POI
     return(out)
}
