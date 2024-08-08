#' Latent moderated strctural equations by Klein and Moosbrugger (2000), the ML approach to nonlinear SEM
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param data_transformations Object containing info on possible data transformations.
#' @param prefix an arbitrary prefix for the data set. This prevents issues when using parallelization and Mplus.
#' @param algorithm algorithm to use. Default to INTEGRATION.
#' @param pathLMS path where (temporal) data and scripts for running LMS using Mplus are stored (using \code{MplusAutomation}). Default to \code{NULL}, then \code{tempdir()} is used.
#' @import stats
#' @importFrom stringr str_split
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @return Returns a \code{data.frame} that includes parameter estimates estimated using LMS.
#' @references Klein, A. G., & Moosbrugger, H. (2000). Maximum likelihood estimation of latent interaction effects with the LMS method. _Psychometrika, 65_(4), 457–474. \doi{10.1007/BF02296338}
#' @export


LMS <- function(lavModel_Analysis, data,
                data_transformations = NULL,
                prefix = 1, pathLMS = tempdir(),
                algorithm = "INTEGRATION")
{
     if(!dir.exists(paste0(pathLMS, "/temp"))) dir.create(paste0(pathLMS, "/temp/"))

     input_file <- paste0(pathLMS, "/temp/", prefix, "_runScript_LMS.inp")
     output_file <- paste0(pathLMS, "/temp/", prefix, "_runScript_LMS.out")
     data_file <- paste0(pathLMS, "/temp/", prefix, "_data_temp.dat")
     lavModel_Analysis_LMS <- lavModel_Analysis

     # Mplus only works with upper cases!
     lavModel_Analysis_LMS$matchLabel <- toupper(lavModel_Analysis_LMS$matchLabel)

     # remove non-fittable coefficients
     lavModel_Analysis_LMS <- lavModel_Analysis_LMS[!((grepl(":", lavModel_Analysis_LMS$lhs) &
                                                            (lavModel_Analysis_LMS$op %in% c("~~", "~1"))) |
                                                     (grepl(":", lavModel_Analysis_LMS$rhs) &
                                                           (lavModel_Analysis_LMS$op %in% c("~~", "~1")))),, drop = FALSE]

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

     write.table(x = data_transformed, file = data_file, col.names = FALSE, row.names = FALSE)
     varnames <- names(data_transformed)

     # measurement model
     temp_measurement <- lavModel_Analysis_LMS[lavModel_Analysis_LMS$op == "=~", c("lhs", "rhs", "start", "fixed"), drop = FALSE]
     if(nrow(temp_measurement)>0)
     {
          Measurement_models <- apply(temp_measurement, 1,
                                      FUN = function(x) paste0(x[1], " BY ", x[2],
                                                               ifelse(test = x[4], yes = "@", no = "*"),
                                                               round(as.numeric(x[3]), 3), ";"))
     }else{
          Measurement_models <- ""
     }


     # nonlinear latent effects
     temp_NL_lv <- lavModel_Analysis_LMS[lavModel_Analysis_LMS$RHSvarType == "funExo", c("rhs")]
     if(length(temp_NL_lv) > 0)
     {
          NL_lv_effects <- unique(temp_NL_lv)
          NL_lv_effects_names <- sapply(stringr::str_split(string = NL_lv_effects, pattern = ":"),
                                        FUN = function(x) paste(paste(x, collapse = "_")))
          NL_lv_effects_XWITH <- sapply(stringr::str_split(string = NL_lv_effects, pattern = ":"),
                                        FUN = function(x) paste(paste(x, collapse = "_"), " | ",
                                                                x[1], " XWITH ", x[2], ";", sep = ""))
     }else{
          NL_lv_effects <- NULL
          NL_lv_effects_names <- ""
          NL_lv_effects_XWITH <- ""
     }


     lavModel_Analysis_LMS_temp <- lavModel_Analysis_LMS
     to_replace <- NL_lv_effects; replacement <- NL_lv_effects_names
     if(length(NL_lv_effects) > 0)
     {
          for(i in seq_along(NL_lv_effects))
          {
               lavModel_Analysis_LMS_temp$lhs[lavModel_Analysis_LMS_temp$lhs == to_replace[i]] <- replacement[i]
               lavModel_Analysis_LMS_temp$rhs[lavModel_Analysis_LMS_temp$rhs == to_replace[i]] <- replacement[i]
          }
     }

     Structural_model <- apply(lavModel_Analysis_LMS_temp[lavModel_Analysis_LMS_temp$op == "~", c("lhs", "rhs", "start", "fixed"),
                                                      drop = FALSE], 1,
                               FUN = function(x) paste0(x[1], " ON ", x[2],
                                                        ifelse(test = x[4], yes = "@", no = "*"),
                                                        round(as.numeric(x[3]), 3), ";"))

     # variances and residual variances
     Variances <- apply(lavModel_Analysis_LMS_temp[lavModel_Analysis_LMS_temp$rhs == lavModel_Analysis_LMS_temp$lhs & lavModel_Analysis_LMS_temp$op == "~~",
                                               c("lhs", "rhs", "start", "fixed", "LHSvarType"), drop = FALSE], 1,
                        FUN = function(x){ out <- paste0(x[1],
                                                         ifelse(test = x[4], yes = "@", no = "*"),
                                                         round(as.numeric(x[3]), 3), ";")
                        if(x[5] == "funExo") return("")
                        out})
     Variances <- Variances[Variances != ""]

     temp <- lavModel_Analysis_LMS_temp[lavModel_Analysis_LMS_temp$rhs != lavModel_Analysis_LMS_temp$lhs & lavModel_Analysis_LMS_temp$op == "~~",
                                    c("lhs", "rhs", "start", "fixed", "LHSvarType", "RHSvarType"), drop = FALSE]
     if(nrow(temp) > 0)
     {
          inds <- apply(temp, 1, function(x) x[5] == "funExo" | x[6] == "funExo")
          temp <- temp[!inds, ,drop = FALSE]
          if(nrow(temp) > 0)
          {
               Covariances <- apply(temp, 1,
                                    FUN = function(x){ out <- paste0(x[1], " WITH ", x[2],
                                                                     ifelse(test = x[4], yes = "@", no = "*"),
                                                                     ifelse(test = is.na(round(as.numeric(x[3]), 3)),
                                                                            yes = "", no =  round(as.numeric(x[3]), 3)), ";")
                                    out})
          }else{
               Covariances <- ""
          }
     }else{
          Covariances <- ""
     }


     # Means and intercepts
     temp_means <- lavModel_Analysis_LMS[lavModel_Analysis_LMS$op == "~1", c("lhs", "rhs", "start", "fixed"), drop = FALSE]
     if(nrow(temp_means)>0)
     {
          Mean_models <- apply(temp_means, 1,
                                      FUN = function(x) paste0("[", x[1],
                                                               ifelse(test = x[4], yes = "@", no = "*"),
                                                               ifelse(test = is.na(round(as.numeric(x[3]), 3)),
                                                                      yes = "", no =  round(as.numeric(x[3]), 3)), "];"))
     }else{
          Measurement_models <- ""
     }


     ### nonlinear effects of ovs missing for now

     LMS_input <- paste0('
TITLE:
     temp', prefix, '

DATA:
     FILE = ', prefix,'_data_temp.dat;

VARIABLE:
     NAMES =\n             ', paste(varnames, sep = "", collapse = "\n             "),';

ANALYSIS:
     TYPE = RANDOM;
     process = 1;
     ALGORITHM = ',algorithm,';

MODEL:
     !Measurement Models
             ', paste(Measurement_models, sep = "", collapse ="\n             "),
                         '\n\n     !Monlinear Effects
             ', paste(NL_lv_effects_XWITH, sep = "", collapse ="\n             "),
                         '\n\n     !Structural Model
             ', paste(Structural_model, sep = "", collapse ="\n             "),
                         '\n\n     !Variances
             ', paste(Variances, sep = "", collapse ="\n             "),
                         '\n\n     !Covariances
             ', paste(Covariances, sep = "", collapse ="\n             "),
                         '\n\n     !Means and Intercepts
             ', paste(Mean_models, sep = "", collapse ="\n             "))

     cat(LMS_input, file = input_file, append = FALSE)
     MplusAutomation::runModels(target = input_file)
     fitMplus <- MplusAutomation::readModels(target = output_file)

     # delete temp files
     file.remove(data_file); file.remove(input_file); file.remove(output_file)

     # check convergence of Mplus
     if(length(fitMplus$warnings) > 0L | length(fitMplus$errors) > 0L){
          stop("LMS had estimation problems.")
     }

     Parameters <- data.frame(fitMplus$parameters$unstandardized)

     # further check convergence and plausibility
     if(is.null(Parameters$se))
     {
          Parameters$est <- NA; Parameters$se <- NA; Parameters$est_se <- NA; Parameters$pval <- NA
     }else{
          if(any(Parameters$se > 20*max(as.numeric(abs(lavModel_Analysis_LMS$start)), na.rm = TRUE)))
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

     return(lavModel_Analysis_LMS)
}

