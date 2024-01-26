#' @import stats

sim_and_fit <- function(n, POI, method,
                        alpha,
                        lavModel, lavModel_Analysis,
                        lavModel_attributes, matrices, data_transformations,
                        prefix,
                        sim_seed,
                        FSmethod = "SL",
                        matchPI =TRUE,
                        PIcentering = "doubleMC",
                        ...){
     set.seed(sim_seed); df_POI <- data.frame("matchLabel" = POI)
     data <- simulateNLSEM(n = n, lavModel = lavModel,
                           lavModel_attributes = lavModel_attributes,
                           matrices = matrices)
     if(tolower(method) == "lms")
     {
          # Mplus only works with upper cases: overwrite df_POI
          df_POI <- data.frame("matchLabel" = toupper(POI))
          fit <- try(LMS(lavModel_Analysis = lavModel_Analysis, data = data,
                         data_transformations = data_transformations, prefix = prefix),
                     silent = TRUE)
     }else if(tolower(method) %in% c("path", "regression", "scaleregression", "sr", "reg")){
          fit <- try(SR(lavModel_Analysis = lavModel_Analysis, data = data,
                        data_transformations = data_transformations), silent = TRUE)
     }else if(tolower(method) %in% c("fsr", "factorscores", "fs")){
          fit <- try(FSR(lavModel_Analysis = lavModel_Analysis, data = data,
                        data_transformations = data_transformations,
                        FSmethod = FSmethod), silent = TRUE)
     }else if(tolower(method) %in% c("upi", "pi", "prodcuctindicator")){
          fit <- try(UPI(lavModel_Analysis = lavModel_Analysis, data = data,
                         data_transformations = data_transformations),
                     silent = TRUE)
     }
     if(!inherits(fit, "try-error"))
     {
          fit_temp <- fit[fit$matchLabel %in% df_POI$matchLabel,,drop = FALSE]
          # sort by POI
          fit <- merge(df_POI, fit_temp, by.x = "matchLabel", sort = FALSE)

          est <- try(fit$est, silent = TRUE); se <- try(fit$se, silent = TRUE)

          if(!inherits(est/se, "try-error"))
          {
               fitOK <- TRUE #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Probably add better check ---------------------------------------------
               if(any(is.na(se))) fitOK <- FALSE
               out <- list("est" = est, "se" = se)

               }else{
               out <- lapply(1:2, function(x) return(rep(NA, length(POI))))
               names(out) <- c("est", "se")
               fitOK <- FALSE
          }
     }else{
          out <- lapply(1:2, function(x) return(rep(NA, length(POI))))
          names(out) <- c("est", "se")
          fitOK <- FALSE
     }
     out <- lapply(out, function(x){names(x) <- POI
                                   return(x)})
     out$fitOK <- fitOK
     return(out)
}
