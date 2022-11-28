#' Latent moderated strctured equations by Klein and Moosbrugger (2000), the ML approach to nonlinear SEM
#' @param lavModel_Analysis the lavModel_Analysis object
#' @param data set to fit
#' @param prefix an arbitrary prefix for the data set. This prevents issues when using parallelization and Mplus.
#' @param algorithm algorithm to use. Default to INTEGRATION.
#' @import MplusAutomation
#' @export




LMS <- function(lavModel_Analysis, data, prefix = 1, algorithm = "INTEGRATION")
{
     if(!dir.exists("temp")) dir.create("temp/")

     input_file <- paste0("temp/", prefix, "_runScript_LMS.inp")
     output_file <- paste0("temp/", prefix, "_runScript_LMS.out")
     data_file <- paste0("temp/", prefix, "data_temp.dat")


     write.table(x = data, file = data_file, col.names = F, row.names = F)
     varnames <- names(data)

     # measurement model
     temp_measurement <- lavModel_Analysis[lavModel_Analysis$op == "=~", c("lhs", "rhs", "start", "fixed"), drop = F]
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
     temp_NL_lv <- lavModel_Analysis[lavModel_Analysis$LHSvarType == "funExo", c("lhs")]
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


     lavModel_Analysis_temp <- lavModel_Analysis
     to_replace <- NL_lv_effects; replacement <- NL_lv_effects_names
     if(length(NL_lv_effects) > 0)
     {
          for(i in seq_along(NL_lv_effects))
          {
               lavModel_Analysis_temp$lhs[lavModel_Analysis_temp$lhs == to_replace[i]] <- replacement[i]
               lavModel_Analysis_temp$rhs[lavModel_Analysis_temp$rhs == to_replace[i]] <- replacement[i]
          }
     }

     Structural_model <- apply(lavModel_Analysis_temp[lavModel_Analysis_temp$op == "~", c("lhs", "rhs", "start", "fixed"),
                                                      drop = F], 1,
                               FUN = function(x) paste0(x[1], " ON ", x[2],
                                                        ifelse(test = x[4], yes = "@", no = "*"),
                                                        round(as.numeric(x[3]), 3), ";"))

     # variances and residual variances
     Variances <- apply(lavModel_Analysis_temp[lavModel_Analysis_temp$rhs == lavModel_Analysis_temp$lhs & lavModel_Analysis_temp$op == "~~",
                                               c("lhs", "rhs", "start", "fixed", "LHSvarType"), drop = F], 1,
                        FUN = function(x){ out <- paste0(x[1],
                                                         ifelse(test = x[4], yes = "@", no = "*"),
                                                         round(as.numeric(x[3]), 3), ";")
                        if(x[5] == "funExo") return("")
                        out})
     Variances <- Variances[Variances != ""]

     temp <- lavModel_Analysis_temp[lavModel_Analysis_temp$rhs != lavModel_Analysis_temp$lhs & lavModel_Analysis_temp$op == "~~",
                                    c("lhs", "rhs", "start", "fixed", "LHSvarType", "RHSvarType"), drop = F]
     if(nrow(temp) > 0)
     {
          inds <- apply(temp, 1, function(x) x[5] == "funExo" | x[6] == "funExo")
          temp <- temp[!inds, ,drop = F]
          if(nrow(temp) > 0)
          {
               Covariances <- apply(temp, 1,
                                    FUN = function(x){ out <- paste0(x[1], " WITH ", x[2],
                                                                     ifelse(test = x[4], yes = "@", no = "*"),
                                                                     round(as.numeric(x[3]), 3), ";")
                                    out})
          }else{
               Covariances <- ""
          }
     }else{
          Covariances <- ""
     }




     ### nonlinear effects of ovs missing for now

     LMS_input <- paste0('
TITLE:
     temp', prefix, '

DATA:
     FILE = ', prefix,'data_temp.dat;

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
             ', paste(Covariances, sep = "", collapse ="\n             "))

     cat(LMS_input, file = input_file, append = F)
     MplusAutomation::runModels(target = input_file)
     fitMplus <- MplusAutomation::readModels(target = output_file)

     fitMplus$parameters$unstandardized

     # delete temp files
     file.remove(data_file); file.remove(input_file); file.remove(output_file)
}

