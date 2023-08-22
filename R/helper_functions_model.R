#' @importFrom stringr str_count
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @import stats

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

     # check variable names
     max_string_length <- max(stringr::str_count(unique(c(lavModel$rhs,
                                                          lavModel$lhs)[!grepl(":", c(lavModel$rhs,
                                                                                      lavModel$lhs))])))
     if(max_string_length>3) warning(paste0("Consider shorter names, especially when using Mplus to fit models.\n",
                                            "Longest name length = ", max_string_length, "."))

     if(any(grepl(pattern = "_l", unique(c(lavModel$rhs,
                                          lavModel$lhs))))) stop("Do not use '_l' within variable names as this is reserved\nfor manifest variables treated as latent variables.")


     # rename higher order terms (quadratic ones)
     lavModel_attributes <- lavaan::lav_partable_attributes(lavModel)
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

          # add fixes (first indicator for now)
          LV <- names(lavModel_attributes$vnames$lv.marker[[1]])
          OV <- lavModel_attributes$vnames$lv.marker[[1]]
          LV <- LV[OV != ""]; OV <- OV[OV != ""]
          lavModel$fixed <- FALSE
          lavModel$fixed[lavModel$lhs %in% LV & lavModel$op == "=~" & lavModel$rhs %in% OV] <- TRUE



          # check lower terms -----
          added_model_syntax <- ""
          for(mod in funRHS)
          {
               lavModel_temp <- lavModel[lavModel$op == "~" & lavModel$rhs == mod,,drop=FALSE]
               for(dv in unique(lavModel_temp$lhs))
               {
                    iv_mods <- stringr::str_split(string = mod, pattern = ":") |> unlist()
                    iv_temp <- lavModel[lavModel$op == "~" & lavModel$lhs == dv,,drop =FALSE]$rhs
                    if(!iv_mods[1] %in% iv_temp)
                    {
                         added_model_syntax <- paste0(added_model_syntax, "\n", paste0(dv,"~0*",iv_mods[1]))
                    }
                    if(iv_mods[1] != iv_mods[2])
                    {
                         if(!iv_mods[2] %in% iv_temp)
                         {
                              added_model_syntax <- paste0(added_model_syntax, "\n", paste0(dv,"~0*",iv_mods[2]))
                         }
                    }
               }
          }

          return(list("lavModel" = lavModel, "added_model_syntax" = added_model_syntax))
     }else{
          stop(paste0("Functional of variables used several times.\nIncluding: ",
                      paste(unique(c(funLHS[problematic_inds_LHS[complete.cases(problematic_inds_LHS)]],
                                     funRHS[problematic_inds_RHS[complete.cases(problematic_inds_RHS)]])),
                            sep = "", collapse = ", "), "."))
     }
}


# compute the matrices necessary to simulate the data

get_matrices <- function(lavModel, lavModel_attributes){

     # if(any(lavModel$RHSvarType == "obs" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     # if(any(lavModel$LHSvarType == "obsEndo" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     if(any(lavModel$ustart[lavModel$op=="~1"][complete.cases(lavModel$ustart[lavModel$op=="~1"])] != 0)){
          stop("Changing the means of latent variables may influence the parameters\nand changes interpretation. This is not recommended.\nMeans for manifest variables should not influence powers,\nhence, this has not been implemented, yet!")
     }
     # ov and lv variables ----
     ov <- lavModel_attributes$vnames$ov[[1]]
     ov <- ov[!grepl(":", ov)]
     ov.nl <- lavModel_attributes$vnames$ov.interaction[[1]]
     ov.y <- lavModel_attributes$vnames$ov.y[[1]]
     lv <- lavModel_attributes$vnames$lv[[1]]
     lv <- lv[!grepl(":", lv)]
     lv.x <- lavModel_attributes$vnames$lv.x[[1]]
     lv.nox <- lavModel_attributes$vnames$lv.nox[[1]]
     lv.y <- lavModel_attributes$vnames$lv.y[[1]]
     lv.nl <- lavModel_attributes$vnames$lv.interaction[[1]]

     eqs.x <- lavModel_attributes$vnames$eqs.x[[1]]
     eqs.y <- lavModel_attributes$vnames$eqs.y[[1]]

     vs <- unique(c(lv.x, lv.nox, lv.y, lv.nl, eqs.x, eqs.y))

     lv.iv <- lv.x[!(lv.x %in% c(lv.nl, lv.nox, eqs.y))]
     ov.iv <- ov[!(ov %in% lavModel$rhs[lavModel$op == "=~"] | ov %in% c(ov.y, eqs.y, ov.nl))]


     # Lambda ----
     if(length(lv) != 0)
     {
          Lambda <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "lambda"]), nrow = max(lavModel$row[lavModel$mat == "lambda"]))
          colnames(Lambda) <- lv
          rownames(Lambda) <- ov[ov %in% lavModel$rhs[lavModel$op == "=~"]] # only measurements that have a measurement model

          RowCol_lambda <- lavModel[lavModel$mat == "lambda",c("row", "col")]
          parameter_starts_lambda <- lavModel$start[lavModel$mat=="lambda"]
          for(i in 1:nrow(RowCol_lambda))
          {
               Lambda[RowCol_lambda[i,1], RowCol_lambda[i,2]] <- parameter_starts_lambda[i]
          }
     }else{
          Lambda <- NULL
     }



     ## theta -----
     if(!is.null(Lambda))
     {
          Theta <- matrix(0, ncol = max(lavModel$col[lavModel$mat == "theta"]), nrow = max(lavModel$row[lavModel$mat == "theta"]))
          RowCol_theta <- lavModel[lavModel$mat == "theta",c("row", "col")]
          colnames(Theta) <- rownames(Theta) <-  rownames(Lambda)
          parameter_starts_theta <- lavModel$start[lavModel$mat=="theta"]
          for(i in 1:nrow(RowCol_theta))
          {
               # symmetric Theta
               Theta[RowCol_theta[i,1], RowCol_theta[i,2]] <- Theta[RowCol_theta[i,2], RowCol_theta[i,1]] <- parameter_starts_theta[i]
          }
     }else{
          Theta <- matrix(0, nrow = length(ov), ncol = length(ov))
          colnames(Theta) <- rownames(Theta) <-  ov
     }


     ## Beta ----
     Beta <- matrix(0, ncol = max(lavModel[lavModel$mat=="beta", c("row", "col")]),
                       nrow = max(lavModel[lavModel$mat=="beta", c("row", "col")]))
     RowCol_beta <- lavModel[lavModel$mat == "beta",c("row", "col")]
     Beta_names <- rep("", max(lavModel[lavModel$mat=="beta", c("row", "col")]))
     Beta_names[lavModel$row[lavModel$mat=="beta" & lavModel$op == "~"]] <- lavModel$lhs[lavModel$mat=="beta" & lavModel$op == "~"]
     Beta_names[lavModel$col[lavModel$mat=="beta" & lavModel$op == "~"]] <- lavModel$rhs[lavModel$mat=="beta" & lavModel$op == "~"]
     Beta_names[lavModel$row[lavModel$mat=="beta" & lavModel$op == "=~"]] <- lavModel$rhs[lavModel$mat=="beta" & lavModel$op == "=~"]
     Beta_names[lavModel$col[lavModel$mat=="beta" & lavModel$op == "=~"]] <- lavModel$lhs[lavModel$mat=="beta" & lavModel$op == "=~"]
     rownames(Beta) <- Beta_names; colnames(Beta) <- Beta_names
     parameter_starts_beta <- lavModel$start[lavModel$mat=="beta"]
     for(i in 1:nrow(RowCol_beta))
     {
          # symmetric Theta
          Beta[RowCol_beta[i,1], RowCol_beta[i,2]]  <- parameter_starts_beta[i]
     }
     # keep rows only if they belong to dependent variables in the model
     Beta_dv <- Beta[!sapply(1:nrow(Beta), function(i) all(Beta[i,]==0)),, drop = FALSE]

     ## Psi ----
     Psi <- matrix(0, ncol = max(lavModel[lavModel$mat == "psi", c("col", "row")]),
                      nrow = max(lavModel[lavModel$mat == "psi", c("col", "row")]))
     RowCol_psi <- lavModel[lavModel$mat == "psi",c("row", "col")]
     Psi_names <- rep("", max(lavModel[lavModel$mat=="psi", c("row", "col")]))
     Psi_names[lavModel$row[lavModel$mat=="psi" & lavModel$op == "~~"]] <- lavModel$lhs[lavModel$mat=="psi" & lavModel$op == "~~"]
     Psi_names[lavModel$col[lavModel$mat=="psi" & lavModel$op == "~~"]] <- lavModel$rhs[lavModel$mat=="psi" & lavModel$op == "~~"]
     rownames(Psi) <- Psi_names; colnames(Psi) <- Psi_names
     parameter_starts_psi <- lavModel$start[lavModel$mat=="psi"]
     for(i in 1:nrow(RowCol_psi))
     {
          # symmetric Psi
          Psi[RowCol_psi[i,1], RowCol_psi[i,2]] <- Psi[RowCol_psi[i,2], RowCol_psi[i,1]] <- parameter_starts_psi[i]
     }
     # remove nl parts (they are implied by lower order moments)
     Psi <- Psi[!grepl(":", Psi_names),  !grepl(":", Psi_names)]

     # full Lambda for simulation ----
     if(length(ov.iv) > 0)
     {
          if(!is.null(Lambda))
          {
               Lambda_full <- rbind(Lambda, matrix(0, nrow = length(ov.iv), ncol = ncol(Lambda)))
               Lambda_full <- cbind(Lambda_full, matrix(0, nrow = nrow(Lambda_full),
                                                        ncol = sum(apply(Theta, 1, function(x) all(x==0))) + length(ov.iv)))
               colnames(Lambda_full)[colnames(Lambda_full) == ""] <- c(colnames(Theta)[apply(Theta, 1, function(x) all(x==0))], ov.iv)
               rownames(Lambda_full)[rownames(Lambda_full) == ""] <- ov.iv
          }else{
               Lambda_full <- matrix(0, nrow = length(ov), ncol = length(ov))
               rownames(Lambda_full) <- colnames(Lambda_full) <- ov
          }
          # fix loading of observed ivs on itself on 1
          pos1 <- unlist(sapply(colnames(Lambda_full), function(x) which(x==rownames(Lambda_full))))
          for(i in 1:length(pos1))
          {
               Lambda_full[pos1[i], names(pos1)[i]] <- 1
          }
     }else{
          Lambda_full <- Lambda
     }

     # compute order ----
     IVs <- c(lv.iv, ov.iv)

     order <- rep(NA, length(vs)); type <- rep(NA, length(vs))
     order[1:length(IVs)] <- 1; type[1:length(IVs)] <- "iv"


     if(is.null(lv.nl))
     {
          templv.nl <- NA
     }else{
          templv.nl <- lv.nl
     }
     if(is.null(ov.nl))
     {
          tempov.nl <- NA
     }else{
          tempov.nl <- ov.nl
     }
     lvov <- IVs; templvov.y <- unique(c(ov.y, lv.y))
     k <- max(order, na.rm = TRUE)+1
     temp_Beta_dv <- Beta_dv
     for(i in 1:nrow(Beta_dv))
     {
          if(length(templv.nl) > 0)
          {
               # check whether nonlinear terms can be formed
               inds <- rowSums(t(sapply(stringr::str_split(string = templv.nl, pattern = ":"),
                                        FUN = function(x) x %in% lvov)), na.rm = TRUE) == 2
               if(any(inds))
               {
                    new <-  templv.nl[inds]
                    lvov <- c(lvov, new)
                    order[is.na(order)][1:length(new)] <- k; k <- k + 1
                    type[is.na(type)][1:length(new)] <- "nl"
                    templv.nl <- templv.nl[!inds]
               }
          }
          if(length(tempov.nl) > 0)
          {
               # check whether nonlinear terms can be formed
               inds <- rowSums(t(sapply(stringr::str_split(string = tempov.nl, pattern = ":"),
                                        FUN = function(x) x %in% lvov)), na.rm = TRUE) == 2
               if(any(inds))
               {
                    new <-  tempov.nl[inds]
                    lvov <- c(lvov, new)
                    order[is.na(order)][1:length(new)] <- k; k <- k + 1
                    type[is.na(type)][1:length(new)] <- "nl"
                    tempov.nl <- tempov.nl[!inds]
               }
          }
          if(nrow(temp_Beta_dv) != 0)
          {
               new <- rownames(temp_Beta_dv)[apply(temp_Beta_dv != 0,
                                                   1, FUN = function(y) all(colnames(temp_Beta_dv[,y]) %in% lvov))]
               if(length(new) == 0){
                    concerningDV <- rownames(temp_Beta_dv)
                    concerningIV <- colnames(temp_Beta_dv[, apply(X = temp_Beta_dv != 0, 2, function(y) any(y))])
                    stop(paste0("Model is circular.\nPlease reconsider your model! You could start with these variable.\nIncludes dependent variables: ",
                                paste(concerningDV, sep = "", collapse = ", "), ".\nIncludes independent variables: ",
                                paste(concerningIV, sep = "", collapse = ", "), ".\n"))
               }
               lvov <- c(lvov, new)
               order[is.na(order)][1:length(new)] <- k; k <- k + 1
               type[is.na(type)][1:length(new)] <- "dv"
               temp_Beta_dv <- temp_Beta_dv[!(rownames(temp_Beta_dv) %in% new),,drop = FALSE]
          }
     }
     Order <- data.frame(order, type, lvov)

     vnames <- list("ov" = ov,
                    "ov.iv" = ov.iv,
                    "ov.nl" = ov.nl,
                    "ov.y" = ov.y,
                    "ov.x" = lavModel_attributes$vnames$ov.x[[1]],
                    "ov.nox" = lavModel_attributes$vnames$ov.nox[[1]],
                    "lv" = lv,
                    "lv.iv" = lv.iv,
                    "lv.nl" = lv.nl,
                    "lv.y" = lv.y,
                    "lv.x" = lv.x,
                    "lv.nox" = lv.nox,
                    "eqs.x" = lavModel_attributes$vnames$eqs.x[[1]],
                    "eqs.y" = lavModel_attributes$vnames$eqs.y[[1]])

     HighestOrders <- compute_highest_order(Beta_dv = Beta_dv, Order = Order)

     return(list("Lambda" = Lambda_full, "Theta" = Theta, "Beta_dv" = Beta_dv, "Psi" = Psi,
                 "vnames" = vnames, "Order" = Order, "HighestOrders" = HighestOrders,
                 "Lambda_small" = Lambda))
}

getModel <- function(lavModel)
{
     temp <- lavModel[!((lavModel$LHSvarType == "funExo" & lavModel$RHSvarType == "funExo" & lavModel$op == "~~") |
                             (lavModel$LHSvarType == "funExo" & lavModel$op == "~1")), c("lhs", "op", "rhs")]
     analysisModel <- paste(apply(temp, 1,
                         FUN = function(string) paste(string, sep = "", collapse = " ")),
                         sep = "", collapse = "\n")
     return(analysisModel)
}


# replace variable names in interactions

replace_variable <- function(var_to_replace, var_to_replace_dv, dv2)
{
     dv2_new <- unique(c(sapply(seq_along(var_to_replace_dv),
                                function(replace_ind) stringr::str_replace_all(string = dv2,
                                                                               pattern = var_to_replace,
                                                                               replacement = var_to_replace_dv[replace_ind]))))
     dv2_new_sorted <- unlist(lapply(lapply(stringr::str_split(dv2_new, pattern = ":"), sort),
                                     function(varnames) paste(varnames, sep = "", collapse = ":")))
     return(dv2_new_sorted) # if "var_to_replace" is not included, the unchanged string (sorted) is returned
}


# compute the highest order in the model
# this can be used to predict run time and to justify certain methods

compute_highest_order <- function(Beta_dv, Order)
{
     VarList <- apply(Beta_dv != 0, 1, function(y) colnames(Beta_dv)[y])
     NL_List <- lapply(VarList, FUN = function(l){l[grepl(pattern = ":", x = l)]})
     included_vars_NL <- lapply(NL_List, FUN = function(nls){unique(unlist(stringr::str_split(nls, ":")))})

     VarList_subs <- list(sort(VarList[[1]]))
     if(length(VarList)>1)
     {
          varnamesY <- names(VarList)
          for(i in 2:length(VarList))
          {
               VarList_subs[[i]] <- VarList[[i]]
               tempnames <- varnamesY[1:(i-1)]
               for(tname in tempnames)
               {
                    VarList_subs[[i]] <- replace_variable(var_to_replace = tname,
                                                          var_to_replace_dv = VarList_subs[names(VarList) == tname][[1]],
                                                          dv2 = VarList_subs[[i]]) |> sort()
               }
          }
          names(VarList_subs) <- names(VarList)
     }

     Orders <- lapply(VarList_subs, FUN = function(l) sapply(stringr::str_split(string = l, pattern = ":"), function(x) length(x)))
     highestOrders <- lapply(Orders, max)
     highestOrder <- max(unlist(highestOrders))
     highestOrder_terms <- unname(unlist(VarList_subs)[which(unlist(Orders) == highestOrder)])

     return(list("highestOrder" = highestOrder, "highestOrder_terms" = highestOrder_terms,
                 "Orders" = Orders, "highestOrders" = highestOrders))
}

# add manifests as latent variables via fixing factor loading to 1 and residual variance to 0.001

add_manifests_as_latent <- function(manifest_po, lavModel_Analysis)
{
     for(man in manifest_po)
     {
          # replace manifest by latent
          lavModel_Analysis$lhs <- stringr::str_replace_all(string = lavModel_Analysis$lhs, pattern = man, replacement = paste0(man, "_l"))
          lavModel_Analysis$rhs <- stringr::str_replace_all(string = lavModel_Analysis$rhs, pattern = man, replacement = paste0(man, "_l"))

     }
     temp_model <- paste(sapply(manifest_po, function(man) paste0(man, "_l =~ ", "1*",man, "\n", man, "~~0.001*",man)),
                         sep = "", collapse = "\n")
     temp_lavModel <- lavaan::lavMatrixRepresentation(lavaan::lavaanify(model = temp_model, meanstructure = TRUE,auto.var = TRUE,
                                                                         auto.cov.lv.x = TRUE, auto.cov.y = TRUE,
                                                                         as.data.frame. = TRUE))
     temp_lavModel <- add_varType(temp_lavModel)
     temp_lavModel <- temp_lavModel[temp_lavModel$op != "~1" & !is.na(temp_lavModel$ustart),, drop = FALSE]
     temp_lavModel$row <- temp_lavModel$col <- temp_lavModel$plabel <- NA
     temp_lavModel$fixed <- TRUE
     temp_lavModel$id <- (nrow(lavModel_Analysis)+1):(nrow(lavModel_Analysis)+nrow(temp_lavModel))
     lavModel_Analysis <- rbind(lavModel_Analysis, temp_lavModel)
     return(lavModel_Analysis)
}


# add covariances among IVs for model estimation

add_covariances_to_lavModel <- function(lavModel_Analysis)
{
     # add covariances among iv
     lavModel_Analysis_attributes <- lavaan::lav_partable_attributes(lavModel_Analysis)
     # ov and lv variables ----
     ov <- lavModel_Analysis_attributes$vnames$ov[[1]]
     ov <- ov[!grepl(":", ov)]
     ov.nl <- lavModel_Analysis_attributes$vnames$ov.interaction[[1]]
     ov.y <- lavModel_Analysis_attributes$vnames$ov.y[[1]]
     lv <- lavModel_Analysis_attributes$vnames$lv[[1]]
     lv <- lv[!grepl(":", lv)]
     lv.x <- lavModel_Analysis_attributes$vnames$lv.x[[1]]
     lv.nox <- lavModel_Analysis_attributes$vnames$lv.nox[[1]]
     lv.y <- lavModel_Analysis_attributes$vnames$lv.y[[1]]
     lv.nl <- lavModel_Analysis_attributes$vnames$lv.interaction[[1]]
     eqs.x <- lavModel_Analysis_attributes$vnames$eqs.x[[1]]
     eqs.y <- lavModel_Analysis_attributes$vnames$eqs.y[[1]]

     lv.iv <- lv.x[!(lv.x %in% c(lv.nl, lv.nox, eqs.y))]
     ov.iv <- ov[!(ov %in% lavModel_Analysis$rhs[lavModel_Analysis$op == "=~"] | ov %in% c(ov.y, eqs.y, ov.nl))]

     IVs <- c(lv.iv, ov.iv)
     included_covs_iv <- lavModel_Analysis[lavModel_Analysis$rhs %in% IVs &
                                                lavModel_Analysis$lhs %in% IVs &
                                                lavModel_Analysis$lhs != lavModel_Analysis$rhs, c("lhs", "rhs"), drop = FALSE]
     if(length(IVs) >= 2)
     {
          temp <- c()
          for(iv1 in seq_along(IVs[-1]))
          {
               for(iv2 in (iv1+1):length(IVs))
               {
                    temp <- rbind(temp, c(IVs[iv1], IVs[iv2]))
               }
          }; temp <- data.frame(temp); names(temp) <- c("V1", "V2")
          temp <- temp[temp$V1 != temp$V2, ]
          if(nrow(included_covs_iv) > 0)
          {
               ind <- apply(temp, 1, function(x) any(apply(included_covs_iv, MARGIN = 1, function(y) (all(x == y) | all(unlist(c(x[2], x[1])) == y)))))
               temp <- temp[!ind, , drop = FALSE]
          }
          if(nrow(temp) > 0)
          {
               m_temp <- paste(apply(temp, 1, function(x) paste(x, collapse = "~~0*")), collapse = "\n")
               temp_lavModel <- lavaan::lavMatrixRepresentation(lavaan::lavaanify(model = m_temp, meanstructure = TRUE,auto.var = TRUE,
                                                                                   auto.cov.lv.x = TRUE, auto.cov.y = TRUE,
                                                                                   as.data.frame. = TRUE))
               temp_lavModel <- add_varType(temp_lavModel)
               temp_lavModel <- temp_lavModel[temp_lavModel$op == "~~" & temp_lavModel$rhs != temp_lavModel$lhs,, drop = FALSE]
               temp_lavModel$row <- temp_lavModel$col <- temp_lavModel$plabel <- NA
               temp_lavModel$fixed <- FALSE
               temp_lavModel$id <- (nrow(lavModel_Analysis)+1):(nrow(lavModel_Analysis)+nrow(temp_lavModel))
               temp_lavModel$start <- temp_lavModel$ustart <- ""
               lavModel_Analysis <- rbind(lavModel_Analysis, temp_lavModel)
          }
     }
     return(lavModel_Analysis)
}

# prepare data for manifest IVs/DVs and change lavModel if necessary
# write data_transformations data.frame

handle_manifests <- function(lavModel, treat_manifest_as_latent = "all")
{

     obs <- unique(lavModel$rhs[lavModel$RHSvarType == "obs"])
     manifest_po <- obs[which(!(obs %in% unique(lavModel$rhs[lavModel$op == "=~"])))] # manifest predictor or outcome
     lavModel_Analysis <- lavModel; data_transformations <- NULL

     if(length(manifest_po) > 0)
     {
          if(tolower(treat_manifest_as_latent) == "all"){
               lavModel_Analysis <- add_manifests_as_latent(manifest_po = manifest_po, lavModel_Analysis = lavModel_Analysis)
          }else if(tolower(treat_manifest_as_latent) == "nl"){
               inds_nl_man <- sapply(manifest_po, function(man) any(grepl(pattern = man,
                                                                          x = unique(lavModel_Analysis$rhs[grepl(pattern = ":",
                                                                                                                 x = lavModel_Analysis$rhs)]))))
               manifest_po_nl <- manifest_po[inds_nl_man]
               if(length(manifest_po_nl)>0)
               {
                    lavModel_Analysis <- add_manifests_as_latent(manifest_po = manifest_po_nl,
                                                                 lavModel_Analysis = lavModel_Analysis)
               }
          }else if(tolower(treat_manifest_as_latent) == "ov"){
               # treat ov as ov, as long as no interaction with lv
               inds_nl_man <- sapply(manifest_po, function(man) any(grepl(pattern = man,
                                                                          x = unique(lavModel_Analysis$rhs[grepl(pattern = ":",
                                                                                                                 x = lavModel_Analysis$rhs)]))))
               manifest_po_nl <- manifest_po[inds_nl_man]
               NLs <- unique(lavModel_Analysis$rhs[grepl(pattern = ":",
                                                         x = lavModel_Analysis$rhs)])
               NLs_split_list <- stringr::str_split(string = NLs, pattern = ":")

               # handle variables that need to be treated as ov
               if(length(manifest_po_nl) == 0) manifest_po_nl <- NULL
               # get variable type
               Variable_type_nl <- t(sapply(NLs_split_list,
                                            FUN = function(nls){temp1 <- ifelse(nls[1] %in% manifest_po_nl, "ov", "lv")
                                            temp2 <- ifelse(nls[2] %in% manifest_po_nl, "ov", "lv")
                                            return(c(temp1, temp2))}))
               colnames(Variable_type_nl) <- c(paste0("var", 1:2, "_type"))
               Variables_nl <- matrix(unlist(NLs_split_list), byrow = TRUE, ncol = 2)
               colnames(Variables_nl) <-  c(paste0("var", 1:2))
               # create type of interaction
               NL_type <- apply(Variable_type_nl, 1, function(x) paste(sort(x), collapse = ":"))


               NL_effects <- cbind(data.frame("NLs" = NLs, "NL_type" = NL_type), Variables_nl, Variable_type_nl)

               # check whether any variable needs to be treated as latent and add to model
               manifest_as_lat <-  unique(c(NL_effects[NL_effects$NL_type != "ov:ov",]$var1[NL_effects[NL_effects$NL_type != "ov:ov",]$var1_type == "lv"],
                                            NL_effects[NL_effects$NL_type != "ov:ov",]$var2[NL_effects[NL_effects$NL_type != "ov:ov",]$var2_type == "lv"]))
               # add ov:lv interactions into lavModel_Analysis object by treating ovs as lvs
               if(length(manifest_as_lat)>0)
               {
                    lavModel_Analysis <- add_manifests_as_latent(manifest_po = manifest_as_lat, lavModel_Analysis = lavModel_Analysis)
               }
               # resume with ovs that do not have ov:lv interactions and treat them as ov
               # generate transformation list
               manifest_po_nl <- manifest_po_nl[!(manifest_po_nl %in% manifest_as_lat)]
               if(length(manifest_po_nl) > 0)
               {
                    OV_NL <- stringr::str_split(string = NL_effects[NL_effects$NL_type == "ov:ov",, drop = FALSE]$NLs, pattern = ":")
                    data_transformations <- c()
                    for(i in seq_along(OV_NL))
                    {
                         data_transformations <- rbind(data_transformations, c(paste(OV_NL[[i]], sep = "", collapse = "_"),
                                                                               OV_NL[[i]][1], OV_NL[[i]][2],
                                                                               paste(OV_NL[[i]], sep = "", collapse = ":")))
                    }; data_transformations <- data.frame(data_transformations); names(data_transformations) <- c("newname", "V1", "V2", "oldname")
               }
          }
     }

     if(!is.null(data_transformations)){ # if is not null, lavModel_Analysis needs to be adapted
          for(i in 1:nrow(data_transformations))
          {
               lavModel_Analysis$RHSvarType[lavModel_Analysis$rhs == data_transformations$oldname[i]] <- "obs"
               lavModel_Analysis$LHSvarType[lavModel_Analysis$lhs == data_transformations$oldname[i]] <- "obs"
               lavModel_Analysis$ustart[lavModel_Analysis$lhs == data_transformations$oldname[i] & lavModel_Analysis$op == "~1"] <-
                    lavModel_Analysis$start[lavModel_Analysis$lhs == data_transformations$oldname[i] & lavModel_Analysis$op == "~1"] <- ""

               lavModel_Analysis$rhs[lavModel_Analysis$rhs == data_transformations$oldname[i]] <- data_transformations$newname[i]
               lavModel_Analysis$lhs[lavModel_Analysis$lhs == data_transformations$oldname[i]] <- data_transformations$newname[i]
          }
     }

     # add covariances among IVs
     lavModel_Analysis <- add_covariances_to_lavModel(lavModel_Analysis = lavModel_Analysis)

     # fix mean structure
     if(sum(lavModel_Analysis$op == "~1" & lavModel_Analysis$LHSvarType != "obs")>0L)
     {
          lavModel_Analysis[lavModel_Analysis$op == "~1" & lavModel_Analysis$LHSvarType != "obs",]$fixed <- TRUE
     }

     return(list("lavModel_Analysis" = lavModel_Analysis,
                 "data_transformations" = data_transformations[!duplicated(data_transformations),,drop = FALSE]))
}






