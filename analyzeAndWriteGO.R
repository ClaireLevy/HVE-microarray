source("getAllDavid.R")
source("getAllInnateDB.R")
# p < 0.05 for days 4, 7, and 14
# df should be allD50 or allD500


analyzeAndWriteGO <- function(df, concentration, filename, direction) {
   source('subsetToOverlappingGoTerms.R')
      
   toWrite <- subsetToOverlappingGoTerms(df, direction, FALSE)

   toWrite<-toWrite %>%
      select(Day,Category,Term,AdjustedPValue)
   
   toWrite <- dcast(toWrite, Category + Term~Day,
         value.var = "AdjustedPValue")
   
   # rename columns because dplyr gets mad at columns named numbers
   colnames(toWrite) <- c("Category", "Term", "d4", "d7", "d14")
   
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
     filter(d4 < 0.05, d7 < 0.05, d14 < 0.05) 
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
                 concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}

analyzeAndWriteGOWithDay1 <- function(df, concentration, filename, direction, day1Sig = FALSE) {
   source('subsetToOverlappingGoTerms.R')
   
   toWrite <- subsetToOverlappingGoTerms(df, direction, TRUE)
   
   toWrite<-toWrite %>%
      select(Day,Category,Term,AdjustedPValue)
   
   toWrite<-dcast(toWrite, Category + Term~Day,
      value.var = "AdjustedPValue")
   
   # rename columns because dplyr gets mad at columns named numbers
   colnames(toWrite) <- c("Category", "Term", "d1", "d4", "d7", "d14")
   
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
      filter(d4 < 0.05, d7 < 0.05, d14 < 0.05)   
   if(day1Sig) {
      toWrite <- filter(toWrite, d1 < 0.05)
   } else {
      # do nothing
   }
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
      concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}