analyzeAndWriteDavid <- function(df, concentration, filename, direction) {
   
   castFunction<-function(data){
      dcast(data, Category + Term~Day,
         value.var = "Benjamini")
   }
   
   # a function to write the DAVID analyzed files
   writeDavidAnalyzed <- function(data, concentration, filename) {
      file <- paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",
         concentration, "/", filename, sep = "")
      write.csv(data, file, row.names=FALSE)
   }
   
   # function that takes a data frame with the following columns, in order, 
   # "Category", "Term", "4", "7", "14" 
   # and returns the rows where p-values are < 0.05 on all three days
   allSig <- function(data) {
      # check data for validity
      if (any(colnames(data) != c("Category", "Term", "4", "7", "14"))) {
         cols <- paste(colnames(data), collapse = ", ")
         stop(paste("allSig expects the column names of data to be: ", 
            "'Category, Term, 4, 7, 14', but the column names of data are '",
            cols, "' so this might be a mistake.", sep = ""))
      }
      
      # rename columns because dplyr gets mad at columns named numbers
      colnames(data) <- c("Category", "Term", "d4", "d7", "d14")
      
      # remove rows where one or more p-value is > 0.05
      data %>%
         filter(d4 < 0.05, d7 < 0.05, d14 < 0.05)   
   }
   
   toWrite <- df[df$direction == direction & df$Day != 1, ]
   toWrite<- toWrite %>%
      group_by(Term)%>%
      filter(n()==3)%>% # because there are 4 days
      arrange(Day,Term,PValue)%>%
      ungroup()

   toWrite<-toWrite %>%
      select(Day,Category,Term,Benjamini)
   
   toWrite<-castFunction(toWrite)
   
   toWrite <- allSig(toWrite)
   writeDavidAnalyzed(toWrite, concentration, filename)
   
   
   
}