getDAVID<-function(day, concentration, direction){

setwd(paste("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = ",concentration,
            sep=""))
  
read.csv(file=paste("DAVID",day, concentration,direction,"csv", sep="."),
              header=TRUE)
}