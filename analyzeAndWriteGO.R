source("getAllDavid.R")
source("getAllInnateDB.R")
# p < 0.05 for days 4, 7, and 14
# df should be allD50 or allD500
analyzeAndWriteGO <- function(df, concentration, filename, direction) {
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

###############################################################################
# DAVID
###############################################################################
# Get data
allD50 <- getAllDavid(50)
allD500 <- getAllDavid(500)

# without Day 1
analyzeAndWriteGO(allD50, 50, "overlap.DAVID.not1.50.DOWN.csv", "DOWN")
analyzeAndWriteGO(allD50, 50, "overlap.DAVID.not1.50.UP.csv", "UP")
analyzeAndWriteGO(allD500, 500, "overlap.DAVID.not1.500.UP.csv", "UP")

# with Day 1
analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.DOWN.csv", "DOWN")
analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.UP.csv", "UP", FALSE)
analyzeAndWriteGOWithDay1(allD500, 500, "overlap.DAVID.500.UP.csv", "UP", TRUE)

###############################################################################
# InnateDB
###############################################################################
# Get data
combinedInnate50<-getAllInnateDB(50) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")

combinedInnate500<-getAllInnateDB(500) %>%
   # rename columns to work with function
   mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value.corrected, 
      Benjamini = PValue, Category = "GO_TERM")

# without Day 1
analyzeAndWriteGO(combinedInnate50, 50, "overlap.Innate.not1.50.DOWN.csv", "DOWN")
analyzeAndWriteGO(combinedInnate50, 50, "overlap.Innate.not1.50.UP.csv", "UP")
analyzeAndWriteGO(combinedInnate500, 500, "overlap.Innate.not1.500.DOWN.csv", "DOWN")
analyzeAndWriteGO(combinedInnate500, 500, "overlap.Innate.not1.500.UP.csv", "UP")

# with Day 1
analyzeAndWriteGOWithDay1(combinedInnate50, 50, "overlap.Innate.50.DOWN.csv", "DOWN", FALSE)
analyzeAndWriteGOWithDay1(combinedInnate50, 50, "overlap.Innate.50.UP.csv", "UP", FALSE)
analyzeAndWriteGOWithDay1(combinedInnate500, 500, "overlap.Innate.500.DOWN.csv", "DOWN", TRUE)
analyzeAndWriteGOWithDay1(combinedInnate500, 500, "overlap.Innate.500.UP.csv", "UP", TRUE)
