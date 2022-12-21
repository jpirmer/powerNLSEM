#' Unconstrained Product Indicator approach by Marsh et al. (2004), with extensions by Kelava and Brandt (2009)
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param data_transformations Data transformations
#' @export


UPI <- function(lavModel_Analysis, data,
                data_transformations = NULL)
{
     lavModel_Analysis_UPI <- lavModel_Analysis

     lavModel_attributes

     ## use semTools for product indicators
     ## use MLR

     # transform data
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



     Parameters <- data.frame(fitMplus$parameters$unstandardized)

     # check convergence
     if(is.null(Parameters$se))
     {
          Parameters$est <- NA; Parameters$se <- NA; Parameters$est_se <- NA; Parameters$pval <- NA
     }else{
          if(any(Parameters$se > 20*max(as.numeric(lavModel_Analysis_LMS$start), na.rm = T)))
          {
               Parameters$est <- NA; Parameters$se <- NA; Parameters$est_se <- NA; Parameters$pval <- NA
          }
     }

     header <- Parameters$paramHeader
     param <- Parameters$param
     Label <- header
     Label <- stringr::str_replace(Label, pattern = ".BY", replacement = "=~")
     Label <- stringr::str_replace(Label, pattern = ".ON", replacement = "~")
     Label <- stringr::str_replace(Label, pattern = ".WITH", replacement = "~~")
     Label[header == "Residual.Variances"] <- paste0(param[header == "Residual.Variances"], "~~")
     Label[header == "Variances"] <- paste0(param[header == "Variances"], "~~")
     Label[header == "Intercepts" | header == "Means"] <- param[header == "Intercepts" | header == "Means"]
     param[header == "Intercepts" | header == "Means"] <- "~1"
     Parameters$matchLabel <- toupper(paste0(Label, param))
     Parameters$matchLabel <- stringr::str_replace_all(string = Parameters$matchLabel, pattern = "_", replacement = ":")

     lavModel_Analysis_LMS <- merge(x = lavModel_Analysis_LMS, y = Parameters, by = "matchLabel")
     lavModel_Analysis_LMS <- lavModel_Analysis_LMS[order(lavModel_Analysis_LMS$id),]
     lavModel_Analysis_LMS$paramHeader <- NULL;      lavModel_Analysis_LMS$param <- NULL;

     # delete temp files
     file.remove(data_file); file.remove(input_file); file.remove(output_file)

     return(lavModel_Analysis_LMS)
}

