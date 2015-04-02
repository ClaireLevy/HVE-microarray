source("getAllDavid.R")
source("getAllInnateDB.R")
# p < 0.05 for days 4, 7, and 14
# df should be allD50 or allD500


analyzeAndWriteGO <- function(df, concentration,
                              filename, direction) {
   source('subsetToOverlappingGoTerms.R')
      
   toWrite <- subsetToOverlappingGoTerms(df,direction, FALSE)
   
  
   
  #Make a df with just Term, Day and Count
  toWrite1<-dplyr::select(toWrite,Category,Term,Day,Count,direction)
  
  #Make another df including PValue and the list total.
  toWrite2<-toWrite %>%
      dplyr::select(Term,Day,Pop.Hits,PValue)
  
  #cast the first df to get a column for each dayCount, a row for each
  #term and the gene count as data
  toWrite1<-dcast(toWrite1,Category+Term+direction ~Day,
                  value.var="Count")
  
  # rename columns because dplyr gets mad at columns named numbers
  colnames(toWrite1) <- c("Category","Term","direction" ,
                          "d4Count", "d7Count", "d14Count")
  
  #cast the second df to get a column for each day and rows for
  #each term and the # of genes in that term. values are Pvalues
  toWrite2 <- dcast(toWrite2, Term + Pop.Hits~Day,
         value.var = "PValue")
  
  # rename columns because dplyr gets mad at columns named numbers
 colnames(toWrite2) <- c("Term","GenesInTerm", "d4.PValue", "d7.PValue", "d14.PValue")
   
  #Merge the two dfs by term
  toWrite<-merge(toWrite1, toWrite2, by = "Term")
  
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
     dplyr::filter(d4.PValue < 0.05, d7.PValue < 0.05,
                   d14.PValue < 0.05) 
   
   toWrite <- arrange(toWrite, d4.PValue, d7.PValue, d14.PValue)
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
                 concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}




analyzeAndWriteGOWithDay1 <- function(df, concentration,
                                      filename, direction,
                                      day1Sig = FALSE) {
  
  
  source('subsetToOverlappingGoTerms.R')
   
   toWrite <- subsetToOverlappingGoTerms(df, direction, TRUE)
   
  #Make a df with just Term, Day and Count
  toWrite1<-dplyr::select(toWrite,Category,Term,Day,Count,direction)
  
  #Make another df including PValue and the list total.
  toWrite2<-toWrite %>%
    dplyr::select(Term,Day,Pop.Hits,PValue)
  
  #cast the first df to get a column for each dayCount, a row for each
  #term and the gene count as data
  toWrite1<-dcast(toWrite1,Category+Term+direction ~Day,
                  value.var="Count")
  
  # rename columns because dplyr gets mad at columns named numbers
  colnames(toWrite1) <- c("Category","Term", "direction","d1Count",
                          "d4Count", "d7Count", "d14Count")
  
  #cast the second df to get a column for each day and rows for
  #each term and the # of genes in that term. values are Pvalues
  toWrite2 <- dcast(toWrite2, Term + Pop.Hits~Day,
                    value.var = "PValue")
  
  # rename columns because dplyr gets mad at columns named numbers
  colnames(toWrite2) <- c("Term","GenesInTerm","d1.PValue",
                          "d4.PValue", "d7.PValue", "d14.PValue")
  
  #Merge the two dfs by term
  toWrite<-merge(toWrite1, toWrite2, by = "Term")
  
 
  
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
      dplyr:: filter(d4.PValue < 0.05, d7.PValue < 0.05, d14.PValue < 0.05)   
   if(day1Sig) {
      toWrite <- dplyr::filter(toWrite, d1.PValue < 0.05)
   } else {
      # do nothing
   }

   toWrite <- arrange(toWrite, d4.PValue, d7.PValue, d14.PValue)
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
      concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}

