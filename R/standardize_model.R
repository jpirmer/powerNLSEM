# models should be filled in the way that if x1:x1 has an effect on Y,
#        then also the effect of x1 should be estimated (although the effect might be zero in the population)
#


small_model <- "
# measurement models
X =~ 1*x1 + 0.8*x2 + 0.7*x3
Y =~ 1*y1 + 0.85*y2 + 0.78*y3
Z =~ 1*z1 + 0.9*z2 + 0.6*z3

# structural models
Y ~ .1*X + .1*Z +  .1*X:Z + .1*x1:x1 +  .1*a:a

# covariances
X~~0.5*Z
"

model <- small_model
lavModel <- lavaan:::lavMatrixRepresentation(lavaan::lavaanify(model = model, meanstructure = T,auto.var = T,
                                                               auto.cov.lv.x = T, auto.cov.y = T,
                                                               as.data.frame. = T))
lavModel <- add_varType(lavModel)
checked_lavModel <- check_model(lavModel)
# modify input model
if(checked_lavModel$added_model_syntax != "")
{
     model <- paste0(model, checked_lavModel$added_model_syntax)
     lavModel <- lavaan:::lavMatrixRepresentation(lavaan::lavaanify(model = model, meanstructure = T,auto.var = T,
                                                                    auto.cov.lv.x = T, auto.cov.y = T,
                                                                    as.data.frame. = T))
     lavModel <- add_varType(lavModel)
     checked_lavModel <- check_model(lavModel)
}
lavModel <- checked_lavModel$lavModel


standardize_model <- function(lavModel) {

     # extract attributes and matrices
     lavModel[lavModel$op == "~" & is.na(lavModel$ustart),c("ustart", "start")] <- c(-Inf)
     lavModel_attributes <- lavaan:::lav_partable_attributes(lavModel)
     matrices <- get_matrices(lavModel, lavModel_attributes)

     # check highest order: for now, only quadratic and interaction models of maximum order 2 are implemented!
     if(matrices$HighestOrders$highestOrder > 2) stop("Standardization for the strutural model has only been implemented\nfor terms of the highest order of 2. Please provide all coefficients when using higher order models.")

     mat <- matrices[1:4]
     Order <- matrices$Order
     Psi <- mat$Psi

     # start off with exogeneous "latent" variables
     Psi00 <- Psi[Order$lvov[Order$order == 1],Order$lvov[Order$order == 1]]

     for(ord in 2:max(Order$order))
     {
          # construct nonlinear effects and append
          if(unique(Order$type[Order$order == ord]) == "nl")
          {
               tempOrder <- Order[Order==ord,, drop = F]; NL <- c()
               for(i in 1:nrow(tempOrder))
               {
                    tempVarNames <- stringr::str_split(pattern = ":", tempOrder$lvov[i])[[1]]
                    NL <- cbind(NL, as.matrix(apply(LV[, tempVarNames, drop = F], 1, prod)))
               }
               colnames(NL) <- tempOrder$lvov
               LV <- cbind(LV, NL)
          }
          # construct dependent variable and append
          if(unique(Order$type[Order$order == ord]) == "dv")
          {
               tempOrder <- Order[Order==ord,, drop = F]

               # catch only exisiting variables
               Zeta_temp <- Zeta[, tempOrder$lvov, drop = F]
               Beta_temp <- mat$Beta_dv[tempOrder$lvov, colnames(LV), drop = F]

               # construct DV from existing LV (and ovs) and "residual" Zeta_temp
               DV <- LV %*% t(Beta_temp) + Zeta_temp
               LV <- cbind(LV, DV)
          }
     }

     #  create manifest variables
     Lambda_temp <- mat$Lambda[matrices$vnames$ov, , drop = F]
     OV <- LV[, colnames(Lambda_temp), drop = F] %*% t(Lambda_temp) + Eps[, matrices$vnames$ov, drop = F]
     OV <- as.data.frame(OV)
     if(appendLVs){
          LV <- as.data.frame(LV)
          OV <- cbind(OV, LV[, !(names(LV) %in% names(OV)), drop = F])
     }
     return(OV)
}

