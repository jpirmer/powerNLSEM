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


get_matrices <- function(lavModel, lavModel_attributes){

     # if(any(lavModel$RHSvarType == "obs" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     # if(any(lavModel$LHSvarType == "obsEndo" & lavModel$op == "~")) stop("Manifest variables can only be measurements for now!\nUsing manifest variables as predictors or outcome\nneeds special considerations on the matrices Beta, Psi,\nLambda and Theta and on the equation by equation simulation!")
     if(any(lavModel$ustart[lavModel$op=="~1"] != 0)) stop("Changing the means of latent variables may influence the parameters\nand changes interpretation. This is not recommended.\nMeans for manifest variables should not influence powers,\nhence, this has not been implemented, yet!")

     # ov variables ----
     ov <- lavModel_attributes$vnames$ov[[1]]
     ov <- ov[!grepl(":", ov)]
     ov.nl <- lavModel_attributes$vnames$ov.interaction[[1]]
     ov.y <- lavModel_attributes$vnames$ov.y[[1]]
     lv <- lavModel_attributes$vnames$lv[[1]]
     lv <- lv[!grepl(":", lv)]

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
     Beta_dv <- Beta[!sapply(1:nrow(Beta), function(i) all(Beta[i,]==0)),, drop = F]

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

     pos1 <- unlist(sapply(colnames(Lambda_full), function(x) which(x==rownames(Lambda_full))))
     for(i in 1:length(pos1))
     {
          Lambda_full[pos1[i], names(pos1)[i]] <- 1
     }

     # compute order ----
     lv.x <- lavModel_attributes$vnames$lv.x[[1]]
     lv.nox <- lavModel_attributes$vnames$lv.nox[[1]]
     lv.y <- lavModel_attributes$vnames$lv.y[[1]]
     lv.nl <- lavModel_attributes$vnames$lv.interaction[[1]]

     eqs.x <- lavModel_attributes$vnames$eqs.x[[1]]
     eqs.y <- lavModel_attributes$vnames$eqs.y[[1]]

     vs <- unique(c(lv.x, lv.nox, lv.y, lv.nl, eqs.x, eqs.y))

     lv.iv <- lv.x[!(lv.x %in% lv.nl)]
     ov.iv <- ov[!(ov %in% lavModel$rhs[lavModel$op == "=~"] | ov %in% c(ov.y, eqs.y))]

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
     k <- max(order, na.rm = T)+1
     temp_Beta_dv <- Beta_dv
     for(i in 1:nrow(Beta_dv))
     {
          if(length(templv.nl) > 0)
          {
               # check whether nonlinear terms can be formed
               inds <- rowSums(t(sapply(stringr::str_split(string = templv.nl, pattern = ":"),
                                        FUN = function(x) x %in% lvov)), na.rm = T) == 2
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
                                        FUN = function(x) x %in% lvov)), na.rm = T) == 2
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
               new <- rownames(temp_Beta_dv)[apply(apply(temp_Beta_dv, 1, FUN = function(x) x != 0),
                                                   2, FUN = function(y) all(colnames(temp_Beta_dv[,y]) %in% lvov))]
               if(length(new) == 0){
                    concerningDV <- rownames(temp_Beta_dv)
                    concerningIV <- colnames(temp_Beta_dv[, apply(X = t(apply(temp_Beta_dv, 1, FUN = function(x) x != 0)), 2, function(y) any(y))])
                    stop(paste0("Model is circular.\nPlease reconsider your model! You could start with these variable.\nIncludes dependent variables: ",
                                paste(concerningDV, sep = "", collapse = ", "), ".\nIncludes independent variables: ",
                                paste(concerningIV, sep = "", collapse = ", "), ".\n"))
               }
               lvov <- c(lvov, new)
               order[is.na(order)][1:length(new)] <- k; k <- k + 1
               type[is.na(type)][1:length(new)] <- "dv"
               temp_Beta_dv <- temp_Beta_dv[!(rownames(temp_Beta_dv) %in% new),,drop = F]
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

     return(list("Lambda" = Lambda_full, "Theta" = Theta, "Beta_dv" = Beta_dv, "Psi" = Psi,
                 "vnames" = vnames, "Order" = Order, "Lambda_small" = Lambda))
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
