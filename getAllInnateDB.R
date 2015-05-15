# Reads in all of the files created by InnateDB for the desired concentration
# and returns them as a data frame




getAllInnateDB <- function(concentration) {
  require(dplyr)
  require(stringr)
   folder <- "J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = "
   folder <- paste(folder, concentration, sep = "")
   #make a list of the files   
   files<-list.files(folder, pattern="^InnateDB")
   
   files <- paste(folder, "/", files, sep = "")
  
   #read in all the files in the list. Now have a list of dfs
   data<-lapply(files,read.csv, sep="\t",quote="")
   
   names(data) <- stringr::str_replace(files, pattern = ".csv",
                                       replacement = "")
   
   #make a list for day and direction 
   dayList<-list(1,1,14,14,4,4,7,7)
   directionList<-list("DOWN","UP")
   
  #add columns to each df in a list of dfs
   data<-Map(cbind,data,Day=dayList,direction=directionList)
   
  #rbind all the data frames in the list together
   suppressWarnings(data <- dplyr::rbind_all(data))
  
   data<-dplyr::select(data,-Pathway.p.value..corrected.)
  
   # rename monstrous columns
   colnames(data)[colnames(data) == 
         "Genes..Symbol.IDBG.ID.Ensembl.Entrez.Fold.Change.P.Value."] <- 
      "Genes.Symbol.IDBG.ID.Ensembl.Entrez.Fold.Change.P.Value"
  
   colnames(data)[5:6]<-c("Count","Pop.Hits")
   # create InnateDB columns that match DAVID columns
   data <- data  %>% mutate(Term = paste(Pathway.Id, Pathway.Name, sep = "~"),
      PValue = Pathway.p.value, Category = Source.Name)
   
   # give columns correct types
   data <- data  %>% mutate(
      Pathway.p.value = as.numeric(Pathway.p.value),
      PValue = as.numeric(PValue))
   
   data
}