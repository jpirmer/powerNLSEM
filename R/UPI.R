#' Unconstrained Product Indicator approach by Marsh et al. (2004), with extensions by Kelava and Brandt (2009)
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param data_transformations Data transformations
#' @param matchPI Logical passed to \code{semTools::indProd} in order to compute the product indicators: Specify TRUE to use match-paired approach (Marsh, Wen, & Hau, 2004). If FALSE, the resulting products are all possible products. Default to \code{TRUE}. The observations are matched by order given when specifying the measurement model.
#' @param PIcentering String indicating which method of centering should be used when constructing product indicators. String is converted to the arguments \code{meanC}, \code{doubleMC}, and \code{residualMC}, of the \code{semTools::indProd} function. Default to \code{"doubleMC"} for double mean centering the resulting products (Lin et. al., 2010). Use \code{"meanC"} for mean centering the main effect indicator before making the products or \code{"residualC"} for residual centering the products by the main effect indicators (Little, Bovaird, & Widaman, 2006). \code{"none"} or any other input than the previously described results in no centering (use with caution!).
#' @param liberalInspection Logical whether the inspection of estimation truthworthiness should be very liberal (i.e., allowing for non-positive definite Hessians in standard error estimation or non-positive residual covariance matrices or latent covariance matrices). Default to \code{FALSE}. Being liberal is not adviced and should be checked for a single data set!
#' @return Returns a \code{data.frame} that includes parameter estimates estimated using UPI.
#' @references Kelava, A., & Brandt, H. (2009). Estimation of nonlinear latent structural equation models using the extended unconstrained approach. _Review of Psychology, 16_(2), 123–132.
#' @references Lin, G. C., Wen, Z., Marsh, H. W., & Lin, H. S. (2010). Structural equation models of latent interactions: Clarification of orthogonalizing and double-mean-centering strategies. _Structural Equation Modeling, 17_(3), 374–391. \doi{10.1080/10705511.2010.488999}
#' @references Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions among latent variables. _Structural Equation Modeling, 13_(4), 497–519. \doi{10.1207/s15328007sem1304_1}
#' @references Marsh, H. W., Wen, Z. & Hau, K. T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. _Psychological Methods, 9_(3), 275–300. \doi{10.1037/1082-989X.9.3.275}
#' @references Marsh, H. W.,  Wen, Z., Hau, K. T., Little, T. D., Bovaird, J. A., & Widaman, K. F. (2007). Unconstrained Structural Equation Models of Latent Interactions: Contrasting Residual- and Mean-Centered Approaches. _Structural Equation Modeling: A Multidisciplinary Journal, 14_(4), 570-580. \doi{10.1080/10705510701303921}
#' @import stats
#' @import utils
#' @export


UPI <- function(lavModel_Analysis, data,
                data_transformations = NULL,
                matchPI =TRUE,
                PIcentering = "doubleMC",
                liberalInspection = FALSE)
{
     lavModel_Analysis_UPI <- lavModel_Analysis

     # transform data
     if(!is.null(data_transformations))
     {
          NL_data <- sapply(1:nrow(data_transformations), FUN = function(d){v1 <- data[, data_transformations$V1[d]]
          v2 <- data[, data_transformations$V2[d]]
          return(scale(v1, center = TRUE, scale = FALSE)*scale(v2, center = TRUE, scale = FALSE))
          })
          NL_data <- data.frame(NL_data); names(NL_data) <- data_transformations$newname
          data_transformed <- cbind(data, NL_data)
     }else{
          data_transformed <- data
     }

     ### generate PI and extend model
     temp_NL_lv <- lavModel_Analysis_UPI[lavModel_Analysis_UPI$RHSvarType == "funExo", c("rhs")]
     dataPI <- c()
     if(length(temp_NL_lv) > 0)
     {
          NL_lv_effects <- unique(temp_NL_lv)
          NL_LV <- stringr::str_split(string = NL_lv_effects, pattern = ":")
          ManifestList <- lapply(NL_LV, function(lvs){manifest1 <- lavModel_Analysis_UPI[lavModel_Analysis_UPI$lhs == lvs[1] &
                                                                                              lavModel_Analysis_UPI$op == "=~", ]$rhs
                                                      manifest2 <- lavModel_Analysis_UPI[lavModel_Analysis_UPI$lhs == lvs[2] &
                                                                                              lavModel_Analysis_UPI$op == "=~", ]$rhs
                                                      return(list(manifest1, manifest2))})
          dataPIList <- lapply(ManifestList, function(man){
               data_temp <- semTools::indProd(data = data,
                                              var1 = man[[1]], var2 = man[[2]],
                                              match = matchPI,
                                              meanC = (PIcentering == "meanC"),
                                              doubleMC = (PIcentering == "doubleMC"),
                                              residualC = (PIcentering == "residualC"))[,-c(1:ncol(data))]
               return(list("dataPI" = data_temp, "PIs" = names(data_temp)))})

          dataPI <- do.call(cbind, lapply(dataPIList, "[[", "dataPI"))
          PIList <- sapply(dataPIList, "[[", "PIs", simplify = FALSE)

          # remove duplicates
          PIList <- lapply(PIList, function(pis){tempI <- stringr::str_split(pis, "\\.")
          dupl <- NULL
          for(ii in 1:(length(tempI)-1))
          {
               if(tempI[[ii]][1] == tempI[[ii]][2]) next
               dupl <- c(dupl, ii + which(sapply(tempI[(ii+1):length(tempI)],
                                                 FUN = function(x) all(tempI[[ii]] %in%  x))))
          }
          if(length(dupl) > 0L) return(pis[-dupl])
          return(pis)})

          # remove duplicates in data
          dataPI <- dataPI[, unlist(PIList)]

          # rename latent Products
          LV_PI <- sapply(NL_lv_effects,
                 function(x) stringr::str_replace_all(string = x, pattern = ":",
                                                      replacement = "_"))

          # Add possible covariances among PIs
          model_covResidPI <- ""
          if(length(PIList) > 1L)
          {
               for(i in 1:(length(PIList)-1))
               {
                    for(j in (i+1):length(PIList))
                    {
                         PIsGrid <- expand.grid(PIList[c(i, j)], stringsAsFactors = FALSE)
                         indCor <- apply(PIsGrid, 1, function(pis){pis1 <- stringr::str_split(string = pis[1],
                                                                                              pattern = "\\.")[[1]]
                                                                   pis2 <- stringr::str_split(string = pis[2], pattern = "\\.")[[1]]
                                                                   return(any(sapply(pis1, FUN = function(x) x %in% pis2)))})
                         # model_covResidPI <- paste(model_covResidPI, "\n",
                         #                           apply(PIsGrid[indCor, ], 1, function(x) paste(x, collapse  = " ~~ ")),
                         #                          collapse = "")
                         model_covResidPI <- c(model_covResidPI, apply(PIsGrid[indCor, ], 1, function(x) paste(x, collapse  = " ~~ ")))
                    }
               }
          }
          # remove duplicates in covariances (if there are any)
          if(length(model_covResidPI) > 1L)
          {
               tempI <- stringr::str_split(model_covResidPI, " ~~ ")
               dupl <- NULL
               for(ii in 1:(length(tempI)-1))
               {
                    dupl <- c(dupl, ii + which(sapply(tempI[(ii+1):length(tempI)],
                                                      FUN = function(x) all(tempI[[ii]] %in%  x))))
               }
               if(length(dupl) == 0L){
                    model_covResidPI <- paste(model_covResidPI,
                                              collapse = "\n")
               }else{
                    model_covResidPI <- paste(model_covResidPI[-dupl],
                                              collapse = "\n")
               }
          }

          # generate measurement models for PIs
          model_PI <- ""
          for(i in seq_along(LV_PI))
          {
               model_PI <- paste0(model_PI, "\n", LV_PI[i], " =~ ",
                                  paste(PIList[[i]], collapse = " + "))
          }

          # combine to PI-model
          model_PI <- paste("# PI measurement model", model_PI,
                            "\n\n\n# PI residual covariances",
                            model_covResidPI, collapse = "\n")

          # append data for PIs
          data_transformed <- cbind(data_transformed, dataPI)
     }else{
          model_PI <- ""
     }

     modelUPI <- paste(apply(lavModel_Analysis_UPI[lavModel_Analysis_UPI$op != "~1",
                                                c("lhs", "op", "rhs")], 1,
                          FUN = function(string) paste(string, sep = "", collapse = " ")),
                    sep = "", collapse = "\n")
     modelUPI <- stringr::str_replace_all(modelUPI, pattern = ":", replacement = "_")
     modelUPI <- paste(modelUPI, "\n\n# PI -----\n\n", model_PI, collapse = "\n")

     fitUPI <- suppressWarnings(lavaan::sem(model = modelUPI,
                                            data = data_transformed,
                                            se = "robust"))

     # check convergence and trustworthiness of results
     if(!(fitUPI@optim$converged)) stop("Error: Model did not converge or showed other issues.")
     if(!liberalInspection &
        !lav_object_post_check(fitUPI))  stop("Error: Model did not converge or showed other issues.")


     Parameters <- lavaan::parameterEstimates(fitUPI)
     Parameters <- Parameters[,1:5]
     Parameters$matchLabel <- apply(Parameters[, 1:3], 1, function(x) paste(x, collapse = ""))
     Parameters$matchLabel <- stringr::str_replace_all(string = Parameters$matchLabel, pattern = "_",
                                                       replacement = ":")
     Parameters <- Parameters[,-(1:3)]

     lavModel_Analysis_UPI <- merge(x = lavModel_Analysis_UPI, y = Parameters, by = "matchLabel",
                                   all.x = TRUE, no.dups = FALSE)
     lavModel_Analysis_UPI <- lavModel_Analysis_UPI[order(lavModel_Analysis_UPI$id),]

     return(lavModel_Analysis_UPI)
}


