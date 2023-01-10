#' plot powerNLSEM object
#' @param out object of class powerNLSEM
#' @param min_num_bins minimal number of bins used for aggregating results. Default to 10.
#' @import ggplot2
#' @import dplyr
#' @export

plot.powerNLSEM <- function(out, min_num_bins = 10)
{
     SUMMARY <- out$SigDecisions %>% group_by(Ns) %>% summarise_all(.funs = function(x) mean(x = x, na.rm = T)) %>% ungroup()
     SUMMARY$Ns_count <- as.numeric(table(out$SigDecisions$Ns))
     NAMES <- names(SUMMARY)
     SUMMARY <- data.frame(SUMMARY)
     names(SUMMARY) <- NAMES

     grouping <- rep(0, length(SUMMARY$Ns_count)); i <- 1; ind_min <- 1
     while(i < length(SUMMARY$Ns_count)+1)
     {
          if(suppressWarnings(min(which(cumsum(SUMMARY$Ns_count[grouping == 0]) >= min_num_bins))) == Inf){
               grouping[grouping == 0] <- i - 1
               break
          }
          ind_max <- min(which(cumsum(SUMMARY$Ns_count[grouping == 0]) >= min_num_bins)) + ind_min - 1
          grouping[ind_min:ind_max] <- i
          ind_min <- ind_max + 1; i <- i + 1
     }
     SUMMARY_agg <- c()
     for(i in unique(grouping))
     {
          SUMMARY_agg <- rbind(SUMMARY_agg,
                               colSums(SUMMARY[grouping == i,, drop = F] *
                                            SUMMARY$Ns_count[grouping == i])/
                                    sum(SUMMARY$Ns_count[grouping == i]))
     }; SUMMARY_agg <- data.frame(SUMMARY_agg); names(SUMMARY_agg) <- NAMES
     SUMMARY_agg$Ns_count <- NULL



     df_long <- reshape(SUMMARY_agg, varying = list(names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"]),
                        direction = "long", v.names = c("Power"),
                        times = names(SUMMARY_agg)[names(SUMMARY_agg) != "Ns"], timevar = "Effect")
     gg <- ggplot(data = df_long, aes(Ns, Power, col = Effect))+geom_point(cex = .1)+
          geom_smooth(method = "loess", formula = "y~x")+
          geom_hline(yintercept = out$power)+geom_vline(xintercept = out$N)
     return(gg)
}
