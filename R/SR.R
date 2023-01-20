#' Scale Regression approach
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param data_transformations Data transformations
#' @export


SR <- function(lavModel_Analysis, data,
                data_transformations = NULL)
{
     lavModel_Analysis_SR <- lavModel_Analysis
     lavModel_measurement_SR <- lavaan:::lav_partable_subset_measurement_model(lavModel_Analysis) |> data.frame()
     lavModel_structural_SR <- lavaan:::lav_partable_subset_structural_model(lavModel_Analysis) |> data.frame()
     # treat all latents as ov as they are collapsed to scales
     lavModel_structural_SR$LHSvarType <- "obs"
     lavModel_structural_SR$RHSvarType[lavModel_structural_SR$op != "~1"] <- "obs"
     temp <- handle_manifests(lavModel = lavModel_structural_SR, treat_manifest_as_latent = "ov")
     data_transformations_latent <- temp$data_transformations
     lavModel_structural_SR <- temp$lavModel_Analysis
     if(!is.null(data_transformations))
     {
          data_transformations <- rbind(data_transformations, data_transformations_latent)
          data_transformations <- data_transformations[!duplicated(data_transformations),, drop = F]
     }else{
          data_transformations <- data_transformations_latent
     }


     dataSR <- c()

     # collapse measurement model
     temp_measurement <- lavModel_measurement_SR[lavModel_measurement_SR$op == "=~", c("lhs", "rhs", "start", "fixed"), drop = F]
     if(nrow(temp_measurement)>0)
     {
          for(lhs in unique(temp_measurement$lhs))
          {
               dataSR <- cbind(dataSR,
                               rowMeans(data[, temp_measurement$rhs[temp_measurement$lhs == lhs],
                                             drop = F], na.rm = T))
          }
          dataSR <- data.frame(dataSR)
          names(dataSR) <- unique(temp_measurement$lhs)
     }

     # transform data
     data <- cbind(data, dataSR)
     if(!is.null(data_transformations))
     {
          NL_data <- sapply(1:nrow(data_transformations), FUN = function(d){v1 <- data[, data_transformations$V1[d]]
          v2 <- data[, data_transformations$V2[d]]
          return(scale(v1, center = T, scale = F)*scale(v2, center = T, scale = F))
          })
          NL_data <- data.frame(NL_data); names(NL_data) <- data_transformations$newname
          data_transformed <- cbind(data, NL_data)
     }else{
          data_transformed <- data
     }

     # fit model
     model <- getModel(lavModel_structural_SR)
     fitSR <- suppressWarnings(lavaan::sem(model = model, data = data_transformed, se = "robust"))
     Parameters <- lavaan::parameterEstimates(fitSR)
     Parameters <- Parameters[,1:5]
     Parameters$matchLabel <- apply(Parameters[, 1:3], 1, function(x) paste(x, collapse = ""))
     Parameters$matchLabel <- stringr::str_replace_all(string = Parameters$matchLabel, pattern = "_",
                                                       replacement = ":")
     Parameters <- Parameters[,-(1:3)]

     lavModel_Analysis_SR <- merge(x = lavModel_Analysis_SR, y = Parameters, by = "matchLabel",
                                   all.x = T, no.dups = F)
     lavModel_Analysis_SR <- lavModel_Analysis_SR[order(lavModel_Analysis_SR$id),]

     return(lavModel_Analysis_SR)
}

