# Reads in all of the files created by DAVID for the desired concentration
# and returns them as a data frame
getAllDavid <- function(concentration) {
   folder <- "J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = "
   folder <- paste(folder, concentration, sep = "")
   
   #extract the david data from the folder where I saved it
   #make a list of the files I want
   Dlist<-list.files(folder, pattern = "^DAVID")
   Dlist <- paste(folder, "/", Dlist, sep = "")
   
   #apply read.csv over the list
   data<-lapply(Dlist,read.csv)
   
   #give the list elements names
   names(data)<- stringr::str_replace(Dlist, pattern = ".csv", replacement = "")
   
   #check the order by looking at names(Ddata)
   #make a list for day and direction 
   dayList<-list(1,1,14,14,4,4,7,7)
   directionList<-list("DOWN","UP")
   
   # add columns to the dfs for day and direction
   data<-Map(cbind,data,Day=dayList)
   data<-Map(cbind,data,direction=directionList)
   
   #make the list of dfs into one df
   suppressWarnings(data <- dplyr::rbind_all(data))
   colnames(data)[colnames(data) == "X."] <- "Percentage"
   
   # add column called Adjusted P Value
   data$AdjustedPValue <- data$Benjamini
   
   #add a column with just the go ids
  data<-mutate(data,Pathway.Id=ifelse(str_detect(data$Term,"GO:")==TRUE,
                             substr(data$Term,1,10),NA))
   
  data
}

# Meaning of columns:
#
# http://www.baderlab.org/Software/EnrichmentMap/UserManual
# Category (DAVID category, i.e. Interpro, sp_pir_keywords, ...)
# Term - Gene set name
# Count - number of genes associated with this gene set, i.e. DEG in this gene set
# Percentage (gene associated with this gene set/total number of query genes)
# P-value - modified Fisher Exact P-value
# Genes - the list of genes from your query set that are annotated to this gene set.
# List Total - number of genes in your query list mapped to any gene set in this ontology
# Pop Hits - number of genes annotated to this gene set on the background list
# Pop Total - number of genes on the background list mapped to any gene set in this ontology.
# Fold enrichment
# Bonferroni
# Benjamini
# FDR