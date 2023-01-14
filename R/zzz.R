powerNLSEMStartupMessage <- function()
{
     msg <- c(paste0("This is a beta version. Please report any bugs! powerNLSEM version ",
          packageVersion("powerNLSEM")),
          "\nType 'citation(\"powerNLSEM\")' for citing this R package in publications.")
     return(msg)
}

.onAttach <- function(lib, pkg)
{
     # unlock .mclust variable allowing its modification
     #unlockBinding(".powerNLSEM", asNamespace("powerNLSEM"))
     # startup message
     msg <- powerNLSEMStartupMessage()
     if(!interactive())
          msg[1] <- paste("Package 'powerNLSEM' version", packageVersion("powerNLSEM"))
     packageStartupMessage(msg)
     invisible()
}
