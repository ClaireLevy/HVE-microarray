

# This file contains the code necessary to convert the longForm data frame into
# .txt and .csv files to enter into InnateDB and DAVID. 
# In both cases, they only include data for genes that were differentially 
# regulated on the day/concentration of interest (i.e. they don't include the
# non-differentially regulated genes for those conditions). 

# get longForm data frame
source("DEG-to-long-form.R")

# Path to save to. Note that there is a dose = 50 and a dose = 500 folder, so
# this path ends with "dose = " and the 50 or 500 and the slashes are pasted
# in inside of the functions
path <- "J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\differentiallyExpressedGenes\\dose = "

#a function to extract and write data from specific days and concentrations
extract<-function(day,concentration,direction){
   x<-longForm[longForm$Day == day &
         
         longForm$Concentration == concentration &
         
         longForm$Direction == direction,
      c(3,6,8)]
   
   file <- paste(day,concentration,direction,"csv",sep=".")
   
   # tab-delimited for InnateDB
   write.table(x,file=paste(path, concentration, "\\", file, sep = ""),
      sep="\t",row.names=FALSE,col.names=FALSE)
   
   # csv for DAVID
   write.csv(x[,1], file= paste(path, concentration, "\\", file, sep = ""),
      row.names=FALSE)
}

extract(1, 50, "UP")
extract(1, 50, "DOWN")
extract(4, 50, "UP")
extract(4, 50, "DOWN")
extract(7, 50, "UP")
extract(7, 50, "DOWN")
extract(14, 50, "UP")
extract(14, 50, "DOWN")

extract(1, 500, "UP")
extract(1, 500, "DOWN")
extract(4, 500, "UP")
extract(4, 500, "DOWN")
extract(7, 500, "UP")
extract(7, 500, "DOWN")
extract(14, 500, "UP")
extract(14, 500, "DOWN")

# alternative way of expressing the same thing
for (day in c(1, 4, 7, 14)) {
   for (concentration in c(50, 500)) {
      for (direction in c("UP", "DOWN")) {
         extract(day, concentration, direction)
      }
   }
}