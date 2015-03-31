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

# Write for DAVID without day 1
allD50 <- getAllDavid(50)
analyzeAndWriteDavid(allD50, 50, "overlap.DAVID.not1.50.UP.csv", "UP")
analyzeAndWriteDavidWithDay1(allD50, 50, "overlap.DAVID.50.DOWN.csv", "DOWN")

allD500 <- getAllDavid(500)
analyzeAndWriteDavid(allD500, 500, "overlap.DAVID.not1.500.UP.csv", "UP")

# Write for InnateDB without day 1
combinedInnate50<-getAllInnateDB(50) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")

analyzeAndWriteGO(combinedInnate50, 50, "overlap.Innate.not1.50.UP.csv", "UP")
analyzeAndWriteGO(combinedInnate50, 50, "overlap.Innate.not1.50.DOWN.csv", "DOWN")

combinedInnate500<-getAllInnateDB(500) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")
analyzeAndWriteGO(combinedInnate500, 500, "overlap.Innate.not1.500.UP.csv", "UP")
analyzeAndWriteGO(combinedInnate500, 500, "overlap.Innate.not1.500.DOWN.csv", "DOWN")