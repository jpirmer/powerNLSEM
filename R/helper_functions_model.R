#' @import stringr
#' @import lavaan

# model helper functions ----

## add a label for the rhs and lhs of the parameter table to aid handling ----
add_varType <- function(lavModel)
{
     lavModel$LHSvarType <- NA; lavModel$RHSvarType <- NA
     lhs <- lavModel$lhs; op <- lavModel$op; rhs <- lavModel$rhs
     varsLHS <- unique(lhs); varsRHS <- unique(rhs)

     for(i in seq_along(varsLHS))
     {
          if(any(lhs == varsLHS[i] & op == "=~"))
          {
               if(any(lhs == varsLHS[i] & op == "~"))
               {
                    LHSvarType <- "latEndo"
               }else{
                    LHSvarType <- "latExo"
               }
          }else{
               if(any(lhs == varsLHS[i] & op == "~"))
               {
                    LHSvarType <- "obsEndo"
               }else{
                    LHSvarType <- "obs"
               }
          }
          if(grepl(pattern = ":",  x = varsLHS[i]))
          {
               LHSvarType <- "funExo"
          }
          lavModel$LHSvarType[lhs == varsLHS[i]] <- LHSvarType
     }

     LHSvarType <- lavModel$LHSvarType

     varsRHS <- varsRHS[varsRHS != ""]
     for(i in seq_along(varsRHS))
     {
          lavModel$RHSvarType[rhs == varsRHS[i]] <- unique(LHSvarType[lhs == varsRHS[i]])
     }

     lavModel$RHSvarType[lavModel$op == "~1"] <- "intercept"

     lavModel[lavModel$op == "==", c("LHSvarType", "RHSvarType")] <- "constraint"

     return(lavModel)
}

## check if model has impossible constraints or if the model can be simulated ----
check_model <- function(lavModel)
{
     if(any(lavModel$group > 1)) stop("Only one group is allowed, for now!")

     if(any(lavModel$LHSvarType == "funExo" & lavModel$user == 1)) stop("Functions of variables cannot be constraint by the user.\nE.g, the variance or the mean of X:Z (the interaction of two variables) is determined by the variances and means of the underlying variables.")
     if(any(lavModel$RHSvarType == "funExo" & lavModel$op == "~~" & lavModel$user == 1)) stop("Functions of variables cannot be constraint by the user.\nE.g, the variance or the mean of X:Z (the interaction of two variables) is determined by the variances and means of the underlying variables.")


     # rename higher order terms (quadratic ones)
     lavModel_attributes <- lavaan:::lav_partable_attributes(lavModel)
     if("NA" %in% lavModel_attributes$vnames$eqs.x[[1]]) stop("Please do not use NA as a name for a latent variable.")
     problematic_names <- lavModel_attributes$vnames$eqs.x[[1]][sapply(stringr::str_split(string = lavModel_attributes$vnames$eqs.x[[1]],
                                                                                          pattern = ":"), function(x) any(x == "NA"))]
     if(length(problematic_names) > 0)
     {
          problematic_names_replace <- unname(sapply(problematic_names, FUN = function(x){
               splitted <- stringr::str_split(x, pattern = ":")[[1]]
               splitted <- splitted[splitted != "NA"]
               paste0(splitted, ":", splitted)
          }))
          for(i in seq_along(problematic_names))
          {
               lavModel$lhs[lavModel$lhs == problematic_names[i]] <- problematic_names_replace[i]
               lavModel$rhs[lavModel$rhs == problematic_names[i]] <- problematic_names_replace[i]
          }
     }


     # check higher order terms
     lhs <- lavModel$lhs; op <- lavModel$op; rhs <- lavModel$rhs
     varsLHS <- unique(lhs); varsRHS <- unique(rhs)
     funLHS <- varsLHS[grepl(pattern = ":", x = varsLHS)]
     funRHS <- varsRHS[grepl(pattern = ":", x = varsRHS)]

     problematic_inds_LHS <- sapply(seq_along(funLHS), FUN = function(i){
          vars <- stringr::str_split(string = funLHS[i], pattern = ":")[[1]]
          fun_new <- paste0(vars[2], ":", vars[1])
          if(fun_new %in% funLHS[-i]){
               return(i)
          }else{return(NA)}})
     problematic_inds_RHS <- sapply(seq_along(funRHS), FUN = function(i){
          vars <- stringr::str_split(string = funRHS[i], pattern = ":")[[1]]
          fun_new <- paste0(vars[2], ":", vars[1])
          if(fun_new %in% funRHS[-i]){
               return(i)
          }else{return(NA)}})

     if(all(is.na(problematic_inds_LHS)) & all(is.na(problematic_inds_RHS))){
          return(lavModel)
     }else{
          stop(paste0("Functional of variables used several times.\nIncluding: ",
                      paste(unique(c(funLHS[problematic_inds_LHS[complete.cases(problematic_inds_LHS)]],
                                     funRHS[problematic_inds_RHS[complete.cases(problematic_inds_RHS)]])),
                            sep = "", collapse = ", "), "."))
     }
}


get_matrices <- function(lavModel){

     if(any(lavModel$RHSvarType == "obs" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     if(any(lavModel$LHSvarType == "obsEndo" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     if(any(lavModel$ustart[lavModel$op=="~1"] != 0)) stop("Changing the means of latent variables may influence the parameters\nand changes interpretation. This is not recommended.\nMeans for manifest variables should not influence powers,\nhence, this has not been implemented, yet!")

     Lambda <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "lambda"]), nrow = max(lavModel$row[lavModel$mat == "lambda"]))
     colnames(Lambda) <- unique(lavModel$lhs[lavModel$mat == "lambda"])
     rownames(Lambda) <- unique(lavModel$rhs[lavModel$mat == "lambda"])

     RowCol_lambda <- lavModel[lavModel$mat == "lambda",c("row", "col")]
     parameter_starts_lambda <- lavModel$start[lavModel$mat=="lambda"]
     for(i in 1:nrow(RowCol))
     {
          Lambda[RowCol_lambda[i,1], RowCol_lambda[i,2]] <- parameter_starts_lambda[i]
     }

     Theta <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "theta"]), nrow = max(lavModel$row[lavModel$mat == "theta"]))
     RowCol_theta <- lavModel[lavModel$mat == "theta",c("row", "col")]
     #names_Theta <- lavModel[lavModel$mat == "theta",c("lhs", "rhs")]
     #rownames(Lambda)[RowCol_theta[apply(names_Theta, 1, function(x) x[1]==x[2]), 1]]
     colnames(Theta) <- rownames(Lambda)
     rownames(Theta) <- rownames(Lambda)
     parameter_starts_theta <- lavModel$start[lavModel$mat=="theta"]
     for(i in 1:nrow(RowCol_theta))
     {
          # symmetric Theta
          Theta[RowCol_theta[i,1], RowCol_theta[i,2]] <- Theta[RowCol_theta[i,2], RowCol_theta[i,1]] <- parameter_starts_theta[i]
     }



     Beta <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "beta"]), nrow = max(lavModel$row[lavModel$mat == "beta"]))
     RowCol_beta <- lavModel[lavModel$mat == "beta",c("row", "col")]
     lv_names <- lavModel_attributes$vnames$lv[[1]]
     rownames(Beta) <- lv_names[1:nrow(Beta)]; colnames(Beta) <- lv_names[1:ncol(Beta)]
     parameter_starts_beta <- lavModel$start[lavModel$mat=="beta"]
     for(i in 1:nrow(RowCol_beta))
     {
          # symmetric Theta
          Beta[RowCol_beta[i,1], RowCol_beta[i,2]]  <- parameter_starts_beta[i]
     }
     # keep rows only if they belong to dependent variables in the model
     Beta_dv <- Beta[!sapply(1:nrow(Beta), function(i) all(Beta[i,]==0)),]


     Psi <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "psi"]), nrow = max(lavModel$row[lavModel$mat == "psi"]))
     RowCol_psi <- lavModel[lavModel$mat == "psi",c("row", "col")]
     colnames(Psi) <- lv_names
     rownames(Psi) <- lv_names
     parameter_starts_psi <- lavModel$start[lavModel$mat=="psi"]
     for(i in 1:nrow(RowCol_psi))
     {
          # symmetric Psi
          Psi[RowCol_psi[i,1], RowCol_psi[i,2]] <- Psi[RowCol_psi[i,2], RowCol_psi[i,1]] <- parameter_starts_psi[i]
     }


     vnames <- list("ov" = rownames(Lambda),
                    "ov.nl" = lavModel_attributes$vnames$ov.interaction[[1]],
                    "ov.y" = lavModel_attributes$vnames$ov.y[[1]],
                    "ov.x" = lavModel_attributes$vnames$ov.x[[1]],
                    "ov.nox" = lavModel_attributes$vnames$ov.nox[[1]],
                    "lv" = colnames(Beta),
                    "lv.nl" = colnames(Beta)[grepl(colnames(Beta), pattern = ":")],
                    "lv.y" = lavModel_attributes$vnames$lv.y[[1]],
                    "lv.x" = lavModel_attributes$vnames$lv.x[[1]],
                    "lv.nox" = lavModel_attributes$vnames$lv.nox[[1]],
                    "eqs.x" = lavModel_attributes$vnames$eqs.x[[1]],
                    "eqs.y" = lavModel_attributes$vnames$eqs.y[[1]])

     return(list("Lambda" = Lambda, "Theta" = Theta, "Beta_dv" = Beta_dv, "Psi" = Psi, "vnames" = vnames))
}
