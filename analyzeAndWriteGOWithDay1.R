# p < 0.05 for days 4, 7, and 14, but also include day 1 - if day1Sig = FALSE, 
# include all rows sig for other days; if day1Sig = TRUE, also require day1 p
# value < 0.05
# df should be allD50 or allD500
analyzeAndWriteGOWithDay1 <- function(df, concentration, filename, direction, day1Sig = FALSE) {
   source('subsetToOverlappingGoTerms.R')
   
   toWrite <- subsetToOverlappingGoTerms(df, direction, TRUE)
   
   toWrite<-toWrite %>%
      select(Day,Category,Term,Benjamini)
   
   toWrite<-dcast(toWrite, Category + Term~Day,
                  value.var = "Benjamini")
 
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

# DAVID with day 1
allD50 <- getAllDavid(50)
analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.UP.csv", "UP", FALSE)
analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.DOWN.csv", "DOWN", FALSE)

allD500 <- getAllDavid(500)
analyzeAndWriteGOWithDay1(allD500, 500, "overlap.DAVID.500.UP.csv", "UP", TRUE)

# InnateDB with day 1
# Write for InnateDB without day 1
combinedInnate50<-getAllInnateDB(50) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")

analyzeAndWriteGOWithDay1(combinedInnate50, 50, "overlap.Innate.50.UP.csv", "UP", FALSE)
analyzeAndWriteGOWithDay1(combinedInnate50, 50, "overlap.Innate.50.DOWN.csv", "DOWN", FALSE)

combinedInnate500<-getAllInnateDB(500) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")
analyzeAndWriteGOWithDay1(combinedInnate500, 500, "overlap.Innate.500.UP.csv", "UP", TRUE)
analyzeAndWriteGOWithDay1(combinedInnate500, 500, "overlap.Innate.500.DOWN.csv", "DOWN", TRUE)
