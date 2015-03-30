# p < 0.05 for days 4, 7, and 14
# df should be allD50 or allD500
analyzeAndWriteDavid <- function(df, concentration, filename, direction) {
   source('subsetToOverlappingGoTerms.R')

      
   toWrite <- subsetToOverlappingGoTerms(df, direction, FALSE)

   toWrite<-toWrite %>%
      select(Day,Category,Term,Benjamini)
   
   toWrite <- dcast(toWrite, Category + Term~Day,
         value.var = "Benjamini")
   
   # rename columns because dplyr gets mad at columns named numbers
   colnames(toWrite) <- c("Category", "Term", "d4", "d7", "d14")
   
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
     filter(d4 < 0.05, d7 < 0.05, d14 < 0.05) 
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
                 concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
   
   
   
}