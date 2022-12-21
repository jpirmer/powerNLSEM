
sim_and_fit <- function(n, POI, method,
                        alpha,
                        lavModel, lavModel_Analysis,
                        lavModel_attributes, matrices, data_transformations,
                        prefix,
                        ...){
     data <- simulateNLSEM(n = n, lavModel = lavModel,
                           lavModel_attributes = lavModel_attributes,
                           matrices = matrices)
     if(tolower(method) == "lms")
     {
          fit <- try(LMS(lavModel_Analysis = lavModel_Analysis, data = data,
                         data_transformations = data_transformations, prefix = prefix), silent = T)
     }
     if(!inherits(fit, "try-error"))
     {
          fit <- fit[fit$matchLabel %in% POI,,drop = F]
          out <- try(abs(fit$est/fit$se), silent = T)
          if(!inherits(out, "try-error"))
          {
               out <- out > qnorm(p = 1-alpha/2)
          }else{
               out <- rep(NA, length(POI))
          }
          names(out) <- POI
     }else{
          out <- rep(NA, length(POI))
     }
     return(out)
}
