#' Factor Score Regression approach
#' @importFrom stringr str_replace_all
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param FSmethod Method to be used to extract factor scores. Default to \code{"SL"} for the Skrondal and Laake approach that uses regression (\code{"regression"}) factor scores for the independendent variables and \code{"Bartlett"} factor scores for the dependent variables.
#' @param data_transformations Data transformations
#' @references **Similar to:** Ng, J. C. K., & Chan, W. (2020). Latent moderation analysis: A factor score approach. _Structural Equation Modeling: A Multidisciplinary Journal, 27_(4), 629–648. <https://doi.org/10.1080/10705511.2019.1664304>.
#' @references Skrondal, A., & Laake, P. (2001). Regression among factor scores. _Psychometrika, 66_(4), 563-575. <https://doi.org/10.1007/BF02296196>
#' @export

FSR <- function(lavModel_Analysis, data, FSmethod = "SL",
                data_transformations = NULL)
{
     lavModel_Analysis_FSR <- lavModel_Analysis
     lavModel_measurement_FSR <- lav_partable_subset_measurement_model(lavModel_Analysis) |> data.frame()
     lavModel_structural_FSR <- lav_partable_subset_structural_model(lavModel_Analysis) |> data.frame()
     # treat all latents as ov as they are collapsed to scales
     lavModel_structural_FSR$LHSvarType <- "obs"
     lavModel_structural_FSR$RHSvarType[lavModel_structural_FSR$op != "~1"] <- "obs"
     temp <- handle_manifests(lavModel = lavModel_structural_FSR, treat_manifest_as_latent = "ov")
     data_transformations_latent <- temp$data_transformations
     lavModel_structural_FSR <- temp$lavModel_Analysis
     if(!is.null(data_transformations))
     {
          data_transformations <- rbind(data_transformations, data_transformations_latent)
          data_transformations <- data_transformations[!duplicated(data_transformations),, drop = FALSE]
     }else{
          data_transformations <- data_transformations_latent
     }

     dataFSR <- c()
     # collapse measurement model
     temp_measurement <- lavModel_measurement_FSR[lavModel_measurement_FSR$op == "=~",
                                                  c("lhs", "rhs", "start", "fixed", "LHSvarType"),
                                                  drop = FALSE]
     if(nrow(temp_measurement)>0)
     {
          for(lhs in unique(temp_measurement$lhs))
          {
               dataFSR <- cbind(dataFSR,
                                getFS(data = data, manifests = temp_measurement$rhs[temp_measurement$lhs == lhs],
                                      latent = lhs, FSmethod = FSmethod,
                                      endogene = unique(temp_measurement$LHSvarType[temp_measurement$lhs == lhs]) == "latEndo"))
          }
          dataFSR <- data.frame(dataFSR)
          names(dataFSR) <- unique(temp_measurement$lhs)

          if(any(apply(dataFSR, 2, function(x) all(is.na(x))))) stop("Error: Factor Scores could not be computed.")
     }

     # transform data
     data <- cbind(data, dataFSR)
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

     # fit model
     model <- getModel(lavModel_structural_FSR)
     fitFSR <- suppressWarnings(lavaan::sem(model = model, data = data_transformed, se = "robust"))
     if(!(lav_object_post_check(fitFSR) & fitFSR@optim$converged)){
          stop("Error: Factor Scores Regression could not be computed.")
     }
     Parameters <- lavaan::parameterEstimates(fitFSR)
     Parameters <- Parameters[,1:5]
     Parameters$matchLabel <- apply(Parameters[, 1:3], 1, function(x) paste(x, collapse = ""))
     Parameters$matchLabel <- stringr::str_replace_all(string = Parameters$matchLabel, pattern = "_",
                                                       replacement = ":")
     Parameters <- Parameters[,-(1:3)]

     lavModel_Analysis_FSR <- merge(x = lavModel_Analysis_FSR, y = Parameters, by = "matchLabel",
                                   all.x = TRUE, no.dups = FALSE)
     lavModel_Analysis_FSR <- lavModel_Analysis_FSR[order(lavModel_Analysis_FSR$id),]

     return(lavModel_Analysis_FSR)
}

