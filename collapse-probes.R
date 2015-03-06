#' Where multiple probe IDs map to the same TargetID, collapse into a single 
#' record per TargetID. E.g. probes 1580204 and 2060433 both detect SFRS18, so 
#' this function collapses the separate entries into one. It depends on the data
#' frame having columns called "TargetID", "ENTREZ_GENE_ID", "Day", 
#' "Concentration", "Log2_fold_change", "Symmetrical_raw_fold_change", 
#' "FDR_adjusted_p_value", and "Direction". 
#' 
#' For the fold changes, it takes either the maximum or the minimum, whichever
#' is bigger in unsigned terms. I.e. if one probe has 1.7 and the other has 
#' -1.8, it will use -1.8. 
#' 
#' For p-value, it takes the minimum.
#' 
#' For direction, it uses UP if the values are all UP or a mixture of UP and 
#' FALSE. Likewise DOWN if all DOWN or a mixture of DOWN and FALSE. It uses
#' FALSE only if all are FALSE. For the rare cases where there are both UP and
#' DOWN it uses a message describing how many probes were up and how many down
#' for that gene.
collapseProbes <- function(dataframe) {
   library(dplyr)
   library(stringr)
   
   # Return the number farther from zero (e.g. c(1.7, -1.8) returns -1.8)
   absMax <- function(data) {
      max <- max(data)
      min <- min(data)
      if (max > abs(min)) {
         return(max)
      } else {
         return(min)
      }
   }
   
   # Return FALSE if all FALSE, UP if all UP or mixture of UP and FALSE, DOWN
   # if all DOWN or mixture of DOWN and FALSE, descriptive message if both UP
   # and DOWN
   getDir <- function(data) {
      if (all(data == "FALSE")) {
         return("FALSE") 
      } else if (all(data %in% c("FALSE", "UP"))) {
         return("UP")
      } else if (all(data %in% c("FALSE", "DOWN"))) {
         return("DOWN")
      } else {
         return(paste("Disagreeing probes:", sum(data == "UP"), "up,", 
            sum(data == "DOWN"), "down"))
      }
   }
      
   dataframe %>%
      group_by(Day, Concentration, TargetID) %>%
      summarize(
         ENTREZ_GENE_ID = ENTREZ_GENE_ID[1],
         Log2_fold_change = absMax(Log2_fold_change),
         Symmetrical_raw_fold_change = absMax(Symmetrical_raw_fold_change),
         FDR_adjusted_p_value = min(FDR_adjusted_p_value),
         Direction = getDir(Direction)
      )
}