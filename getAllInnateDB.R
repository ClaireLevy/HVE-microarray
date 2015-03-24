# Reads in all of the files created by InnateDB for the desired concentration
# and returns them as a data frame
getAllInnateDB <- function(concentration) {

   folder <- "J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = "
   folder <- paste(folder, concentration, sep = "")
      
   files<-list.files(folder, pattern="^InnateDB")
   
   data<-lapply(files,read.csv, sep="\t",quote="")
   
   names(data) <- stringr::str_replace(files, pattern = ".csv", replacement = "")
   
   #make a list for day and direction 
   dayList<-list(1,1,14,14,4,4,7,7)
   directionList<-list("DOWN","UP")
   
   data<-Map(cbind,data,Day=dayList,direction=directionList)
   
   data <- dplyr::rbind_all(data)
   
   # rename monstrous columns
   colnames(data)[colnames(data) == "Pathway.p.value..corrected."] <- 
      "Pathway.p.value.corrected"
   colnames(data)[colnames(data) == 
         "Genes..Symbol.IDBG.ID.Ensembl.Entrez.Fold.Change.P.Value."] <- 
      "Genes.Symbol.IDBG.ID.Ensembl.Entrez.Fold.Change.P.Value"
   
   data
}