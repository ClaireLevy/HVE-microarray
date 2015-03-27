# p < 0.05 for days 4, 7, and 14, but also include day 1 - if day1Sig = FALSE, 
# include all rows sig for other days; if day1Sig = TRUE, also require day1 p
# value < 0.05
# df should be allD50 or allD500
analyzeAndWriteDavidWithDay1 <- function(df, concentration, filename, direction, day1Sig = FALSE) {
   source('subsetToOverlappingGoTerms.R')
   castFunction<-function(data){
      dcast(data, Category + Term~Day,
         value.var = "Benjamini")
   }
   
   # a function to write the DAVID analyzed files
   writeDavidAnalyzed <- function(data, concentration, filename) {
      file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
         concentration, "/", filename, sep = "")
      write.csv(data, file, row.names=FALSE)
   }
   
   # function that takes a data frame with the following columns, in order, 
   # "Category", "Term", "4", "7", "14" 
   # and returns the rows where p-values are < 0.05 on all three days
   allSig <- function(data) {
      # check data for validity
      if (any(colnames(data) != c("Category", "Term", "1", "4", "7", "14"))) {
         cols <- paste(colnames(data), collapse = ", ")
         stop(paste("allSig expects the column names of data to be: ", 
            "'Category, Term, 4, 7, 14', but the column names of data are '",
            cols, "' so this might be a mistake.", sep = ""))
      }
      
      # rename columns because dplyr gets mad at columns named numbers
      colnames(data) <- c("Category", "Term", "d1", "d4", "d7", "d14")
      
      # remove rows where one or more p-value is > 0.05
      data <- data %>%
         filter(d4 < 0.05, d7 < 0.05, d14 < 0.05)   
      if(day1Sig) {
         return(filter(data, d1 < 0.05))
      } else {
         return(data)
      }
   }
   
   toWrite <- subsetToOverlappingGoTerms(df, direction, TRUE)
   
   toWrite<-toWrite %>%
      select(Day,Category,Term,Benjamini)
   
   toWrite<-castFunction(toWrite)
   
   toWrite <- allSig(toWrite)
   writeDavidAnalyzed(toWrite, concentration, filename)
   
   
   
}