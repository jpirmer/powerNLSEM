#' @import stats

sim_and_fit <- function(n, POI, method,
                        alpha,
                        lavModel, lavModel_Analysis,
                        lavModel_attributes, matrices, data_transformations,
                        prefix,
                        sim_seed,
                        FSmethod = "SL",
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
          true <- try(fit$ustart, silent = TRUE)
          pvalue_onesided <- try(pnorm(sign(true)*est / se, lower.tail = FALSE), silent = TRUE)
          pvalue_twosided <- try(2*pnorm(abs(est) / se, lower.tail = FALSE), silent = TRUE)

          if(!inherits(pvalue_onesided, "try-error") & !inherits(pvalue_twosided, "try-error"))
          {
               sigs_onesided <- pvalue_onesided < alpha
               sigs_twosided <- pvalue_twosided < alpha

               fitOK <- TRUE #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Probably add better check ---------------------------------------------

               out <- list("est" = est, "se" = se, "pvalue_onesided" = pvalue_onesided,
                           "pvalue_twosided" = pvalue_twosided,
                           "sigs_onesided" = sigs_onesided, "sigs_twosided" = sigs_twosided)

               }else{
               out <- lapply(1:6, function(x) return(rep(NA, length(POI))))
               names(out) <- c("est", "se", "pvalue_onesided", "pvalue_twosided",
                               "sigs_onesided", "sigs_twosided")
               fitOK <- FALSE
          }
     }else{
          out <- lapply(1:6, function(x) return(rep(NA, length(POI))))
          names(out) <- c("est", "se", "pvalue_onesided", "pvalue_twosided",
                          "sigs_onesided", "sigs_twosided")
          fitOK <- FALSE
     }
     out <- lapply(out, function(x){names(x) <- POI
                                   return(x)})
     out$fitOK <- fitOK
     return(out)
}
