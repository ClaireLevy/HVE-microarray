source("getAllDavid.R")
source("getAllInnateDB.R")
# p < 0.05 for days 4, 7, and 14
# df should be allD50 or allD500



day1Sig<-FALSE
df<-allD50
concentration<-50
filename<-"bk1.csv"
direction<-"UP"

analyzeAndWriteGO <- function(df, concentration,
                              filename, direction) {
   source('subsetToOverlappingGoTerms.R')
      
   toWrite <- subsetToOverlappingGoTerms(df, direction, FALSE)
   
   ######## Finding # of genes assoc with GO ids
   require(topGO)
   
   # I want to make a column in a df of overlap data with the
   #number of genes in each go id.
  
   
   #make a list of unique go Ids from allD50
  List<-as.list(unique(df$Pathway.Id))
   
   #mapping of Go ids to entrez ids for my universe
   #Using topGO functions and methods
   GOID2Gene<-annFUN.org(c("CC","BP","MF"),
                         feasibleGenes = NULL,
                         mapping="org.Hs.eg.db",
                         ID = "entrez")
   
   #geneNames are the GO ids associated with entrez Ids
   #they are the names of the elemnts in the GOID2Gene list
   geneNames<-names(GOID2Gene)
   
   #filter the elements in the universe that overlap
   #with list of GO ids from allD50.
   geneList<-GOID2Gene[geneNames %in% List]
   
   #here is a named list where the names are GO ids from allD50
   #and the element is the number of associated genes
   GOidLength<-lapply(geneList,length)
   
   #convert the GOidLength list into a dataframe
   #used dplyr's as_data_frame so the GO ids aren't turned into
   #rownames
   GOidLength<-as_data_frame(GOidLength)
   
   #now it is wide and I want it to be long so I will melt
   GOidLength<-melt(GOidLength)
   
   #cleanup so it will merge well
   colnames(GOidLength)<-c("Pathway.Id","GenesInGOid")
   GOidLength$Pathway.Id<-as.character(GOidLength$Pathway.Id)
   GOidLength$GenesInGOid<-as.integer(GOidLength$GenesInGOid)
   
   #merge it with the allD50 dataframe
  toWrite<-merge( GOidLength,toWrite, by= "Pathway.Id")
   
  
   ### NOTE: This no longer contains data from non-GO databases
   ### like SPIR and etc.
   
  
  #Make a df with just Term, Day and Count
  toWrite1<-dplyr::select(toWrite,Term,Day,Count)
  
  #Make another df with all columns EXCEPT Count
  toWrite2<-toWrite %>%
      dplyr::select(Day,Category,Term,AdjustedPValue,
                    GenesInGOid)
  #cast the first df to get a column for each dayCount, a row for each
  #term and the gene count as data
  toWrite1<-dcast(toWrite1,Term~Day, value.var="Count")
  
  # rename columns because dplyr gets mad at columns named numbers
  colnames(toWrite1) <- c("Term", "d4Count", "d7Count", "d14Count")
  
  #cast the second df to get a column for each day and rows for
  #each term and the # of gens in that term. values are adjP
  toWrite2 <- dcast(toWrite2, Category + Term + GenesInGOid~Day,
         value.var = "AdjustedPValue")
  
  # rename columns because dplyr gets mad at columns named numbers
 colnames(toWrite2) <- c("Category", "Term","GenesInGOid", "d4", "d7", "d14")
   
  #Merge the two dfs by term
  toWrite<-merge(toWrite1, toWrite2, by = "Term")
  
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
     dplyr::filter(d4 < 0.05, d7 < 0.05, d14 < 0.05) 
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
                 concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}




analyzeAndWriteGOWithDay1 <- function(df, concentration,
                                      filename, direction,
                                      day1Sig = FALSE) {
  
  
  source('subsetToOverlappingGoTerms.R')
   
   toWrite <- subsetToOverlappingGoTerms(df, direction, TRUE)
   
   
   ######## Finding # of genes assoc with GO ids
   require(topGO)
   
   # I want to make a column in a df of overlap data with the
   #number of genes in each go id.
   
   
   #make a list of unique go Ids from allD50
   List<-as.list(unique(df$Pathway.Id))
   
   #mapping of Go ids to entrez ids for my universe
   #Using topGO functions and methods
   GOID2Gene<-annFUN.org(c("CC","BP","MF"),
                         feasibleGenes = NULL,
                         mapping="org.Hs.eg.db",
                         ID = "entrez")
   
   #geneNames are the GO ids associated with entrez Ids
   #they are the names of the elemnts in the GOID2Gene list
   geneNames<-names(GOID2Gene)
   
   #filter the elements in the universe that overlap
   #with list of GO ids from allD50.
   geneList<-GOID2Gene[geneNames %in% List]
   
   #here is a named list where the names are GO ids from allD50
   #and the element is the number of associated genes
   GOidLength<-lapply(geneList,length)
   
   #convert the GOidLength list into a dataframe
   #used dplyr's as_data_frame so the GO ids aren't turned into
   #rownames
   GOidLength<-as_data_frame(GOidLength)
   
   #now it is wide and I want it to be long so I will melt
   GOidLength<-melt(GOidLength)
   
   #cleanup so it will merge well
   colnames(GOidLength)<-c("Pathway.Id","GenesInGOid")
   GOidLength$Pathway.Id<-as.character(GOidLength$Pathway.Id)
   GOidLength$GenesInGOid<-as.numeric(GOidLength$GenesInGOid)
   
   #merge it with the allD50 dataframe
   toWrite<-merge(GOidLength,toWrite, by = "Pathway.Id")
   
   ### NOTE: This no longer contains data from non-GO databases
   ### like SPIR and etc.
   #Make a df with just Term, Day and Count
   toWrite1<-dplyr::select(toWrite,Term,Day,Count)
   
   #Make another df with all columns EXCEPT Count
   toWrite2<-toWrite %>%
     dplyr::select(Day,Category,Term,AdjustedPValue,
                   GenesInGOid)
   #cast the first df to get a column for each dayCount, a row for each
   #term and the gene count as data
   toWrite1<-dcast(toWrite1,Term~Day, value.var="Count")
   
   # rename columns because dplyr gets mad at columns named numbers
   colnames(toWrite1) <- c("Term","d1Count","d4Count", "d7Count", "d14Count")
   
   #cast the second df to get a column for each day and rows for
   #each term and the # of gens in that term. values are adjP
   toWrite2 <- dcast(toWrite2, Category + Term + GenesInGOid~Day,
                     value.var = "AdjustedPValue")
   
   # rename columns because dplyr gets mad at columns named numbers
   colnames(toWrite2) <- c("Category", "Term","GenesInGOid","d1", "d4", "d7", "d14")
   
  #Merge the two dfs by term
  toWrite<-merge(toWrite1, toWrite2, by = "Term")
   
  
  
  
  
   # remove rows where one or more p-value is > 0.05
   toWrite <- toWrite %>%
      dplyr:: filter(d4 < 0.05, d7 < 0.05, d14 < 0.05)   
   if(day1Sig) {
      toWrite <- filter(toWrite, d1 < 0.05)
   } else {
      # do nothing
   }
   
   file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
      concentration, "/", filename, sep = "")
   write.csv(toWrite, file, row.names=FALSE)
}