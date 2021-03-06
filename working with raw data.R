#Working with raw illumina data
#started 26March15


require(dplyr)
library(lumi)
library(limma)

setwd("J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\microarrayData")
#LMF sent FinalReport_exvivoTNF_HVE_032615.txt which is raw
#data with control PROBE profile data rather than control
#GENE profile data, which we tried before and couldn't get lumi
#to read.

#she also sent her code for reading in and normalising (below)

######################## READ IN RAW DATA #########################
#LMF:load data into memory 
RAW<-'FinalReport_032615_exvivo_TNF_HVE.txt'
##LMF:read raw data into lumibatch object
RAW.lumi =lumiR(RAW,
                detectionTh = 0.05, na.rm = TRUE,
                convertNuID = FALSE, dec = '.',
                parseColumnName = FALSE, checkDupId = FALSE,
                QC = TRUE, 
                columnNameGrepPattern = list(exprs='AVG_SIGNAL',
                                             se.exprs='BEAD_STDERR',
                                             detection='Detection Pval',
                                             beadNum='Avg_NBEADS'),
                inputAnnotation=TRUE,
                annotationColumn=c('ILMN_GENE', 'ENTREZ_GENE_ID', 'GI', 'ACCESSION', 'SYMBOL', 'PROBE_ID', 'PROBE_START', 'PROBE_SEQUENCE', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES'),
                verbose = TRUE) 



#all the feature data incl all the various probe ids and symbols
fData<-fData(RAW.lumi)
#just probeID and symbols, remove the rownames, which are the same
#as the ProbeID (different from PROBE_ID)
#use this later for annotation

ProbeIDTargetIDEntrez<-fData%>%
  dplyr::select(ProbeID, TargetID,ENTREZ_GENE_ID)
colnames(ProbeIDTargetIDEntrez)<-c("Probe.ID","TargetID","ENTREZ_GENE_ID")

rownames(ProbeIDTargetIDEntrez)<-NULL

#################### NORMALISATION AND QC ###########################

##LMF:perform vst transformation and rsn normalization
#lumiT returns a LumiBatch object with transformed exprs values
#does log2, vst and cubicRoot transformations
RAWlumi.T = lumiT(RAW.lumi)

#between chip normalization, makes bgrnd same for all?
RAWlumi.N = lumiN(RAWlumi.T, method = "rsn")

#QC of normalized data
RAWlumi.N.Q = lumiQ(RAWlumi.N, detectionTh = 0.05)
summary(RAWlumi.N.Q,"QC")#mean expression is much more uniform

########################## FILTER DETECTION ###############################

#I only want probes where the detection is <0.05 for all 3 donors of dose=0,
#OR all of dose=50 OR all of dose=500 

#from the detection data of the lumiBatch created above:

#For each row (probe) determine whether the detection p values are < 0.05
# for ALL the donors for at least one of the 12 combinations of day+dose

# The resulting logical vector is ALL the probes and the value TRUE or
# FALSE for each of them indicating whether or not  that probe met the condition above



# The resulting logical vector is all the probes and the value TRUE or
# FALSE for each of them telling me if at least 1/12 conditions fits the
#rule that rowSums==3.

allDays.N.Q <- RAWlumi.N.Q

# indices are for HVE_A*, HVE_B*, and HVE_C*, where * is the same number
expressedLogical<-rowSums(detection(allDays.N.Q[,c(2, 18, 14)]) <0.05)==3|
  rowSums(detection(allDays.N.Q[,c(19, 4, 27)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(24, 28, 5)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(31, 22, 32)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(7, 13, 1)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(16, 17, 21)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(35, 9, 29)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(25, 10, 26)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(8, 6, 15)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(33, 11, 12)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(36, 20, 3)]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,c(30, 34, 23)]) <0.05) ==3
  
#Now I need to subset the lumibatch to just contain those TRUE probes.
#To do that, I pass this logical vector to the lumiBatch 
# which subsets the lumiBatch to contain only the probes that are TRUE.
expressed<-allDays.N.Q[expressedLogical,]

dims(RAWlumi.N.Q)
dims(expressed)
#20117 probes left
dims(RAWlumi.N.Q)-dims(expressed)
#removed 27206 probes

############################ FILTER SD  ##################

library(limma)
library(genefilter)

dataMatrix <- exprs(expressed)

#calculate the row standard deviations
rowSdsMatrix<-rowSds(dataMatrix)

#subset dataMatrix to get just the probes where the rowSD is >0.2
dataMatrixSDfilter<-dataMatrix[rowSdsMatrix>=0.2,]

#compare the filtered and unfiltered lumiBatch

dim(dataMatrix)[1]-dim(dataMatrixSDfilter)[1]
#got rid of 9696 probes
dim(dataMatrixSDfilter)[1]
#10421 left

################ FIT MODEL TO ALL DATA GROUPS ################################

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData")

library(limma)
library(genefilter)
library(dplyr)
library(ggplot2)
library(stringr)
load("dataMatrixSDfilter.Rda")# data filtered by SD and detection

#sample key has sample information incl a "Group" column
#which assigns a number to a day + treatment + paired control (8)

# Below is code to make a list of data frames; one for each
#sample group.
# A group is: day + treatment + paired control. Incl all 3 donors.
#Group 1 is day 1, dose = 50 and dose =0
#Group 2 is day 1 dose = 500 and dose = 0
# etc up to group 8 (day 14, dose =500 and dose = 0)


###################### SUBSET LIMMA ################################
load("SampleKey.Rda")

#Function to select a particular group from the limma data set
SelectDataGroup<-function(Group){
  x<-SampleKey[SampleKey$Group == Group,]
  dataMatrixSDfilter[,as.character(x$SampleName)]
}

# I want a list consisting of data frames for each group
#Here is a list of the groups
GroupList<-list(1,2,3,4,5,6,7,8)

#pass the list to the SelectDataGroup function
GroupDataframes<-lapply(GroupList,FUN=SelectDataGroup)


#This file has the code for making the design matrix and targets df
source("DesignMatrixAndTargets.R")


########################### FIT LIMMA ################################

#fit the model using the design matrix and the data frames for each

#group
fit<-lapply(GroupDataframes,FUN=lmFit,design = design)

fit<-lapply(fit, FUN=eBayes)

#generate a top table from the treatment data with BH adj pvalues
#give all entries (so, as many rows as are in fit, number= inf)
#The default topTable only gives the "top" results, but I want 
#all of them.

TT<-lapply(fit,FUN=topTable,coef="Treatmentdrug",
           adjust="BH", number = Inf)



######################## FILTER FOR LOG FC AND PVAL ###############

#filter for logFC
LogFCfiltered<-lapply(TT,subset,logFC>=0.5 | logFC<=-0.5 )

#filter for adj p value
PValLogFCfiltered<-lapply(LogFCfiltered,subset,adj.P.Val <=0.05)


#change the rowsnames(probes) to an actual column
#I think there is a smarter way to do this but this works
#maybe
#lapply(PVallogFCfiltered, rownames,
# PVaLogFCfiltered$Probe.ID)



PValLogFCfiltered<-lapply(GroupList, function(df){
  c<-PValLogFCfiltered[[df]]
  c$Probe.ID<-rownames(c)
  rownames(c)<-NULL
  return(c)
  
})



lapply(PValLogFCfiltered,FUN=nrow)




######################GET HGNC SYMBOLS #################################


setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData")

load("PValLogFCfiltered14Sept15.Rda")



#merge in the symbols names using ProbeIDTargetIDEntrez
#df that I made when I read in the raw data above
load("ProbeIDTargetIDEntrez.Rda")

PValLogFCfiltered<-lapply(PValLogFCfiltered,merge,
                          ProbeIDTargetIDEntrez, by="Probe.ID")

############ COMPARING PROBES UP AND DOWN PER CONCENTRATION#################

#add a column to each df in the list that says which day
#and concentration the data is about

dayList<-list(1,1,4,4,7,7,14,14)

PValLogFCfiltered<-Map(cbind,PValLogFCfiltered,Day=dayList)

concList<-list(50,500,50,500,50,500,50,500)

PValLogFCfiltered<-Map(cbind,PValLogFCfiltered,Concentration=concList)


#unlist the list so there is one big df in longform with the data
#for plotting
#I don't know why this works
processedData<-do.call(rbind, lapply(PValLogFCfiltered,data.frame))

#add a column for the direction of expression (up or down)

processedData<-mutate(processedData,Direction=ifelse(processedData$logFC<0,
                                            "DOWN","UP"))
processedData$Concentration<-as.factor(processedData$Concentration)

processedData$Day<-as.factor(processedData$Day)

probeCountPlot<-ggplot(processedData,aes(x = Day))+
  geom_bar(aes(fill = Direction),position="dodge")+
  scale_fill_manual(values = c("UP"="coral", "DOWN"= "cornflowerblue"))+
  facet_wrap(~Concentration)+
  ggtitle("Count of up and down regulated probes")

ggsave("probeCountPlot.png",dpi=600)

ggplot(processedData,aes(ENTREZ_GENE_ID, Probe.ID ))+
  geom_tile(aes(fill = AveExpr))+
  scale_fill_gradient(low = "white",high = "steelblue")

######################### COLLAPSE VALUES FOR MULTIPLE PROBES #######################

#I am using a slightly modified version of SH's code in "collapse-probes.R" that
#he wrote for use with the CyberT data from LMF. I don't have the 
#"symmetrical raw fold change" column

CLcollapseProbes <- function(dataframe) {
  library(dplyr)
  library(stringr)
  
  # Return the number farther from zero (e.g. c(1.7, -1.8) returns -1.8)
  absMax <- function(data) {
    max <- max(data)
    min <- min(data)
    if (max > abs(min)) {
      return(max)
    } else {
      return(min)
    }
  }
  
  # Return FALSE if all FALSE, UP if all UP or mixture of UP and FALSE, DOWN
  # if all DOWN or mixture of DOWN and FALSE, descriptive message if both UP
  # and DOWN
  getDir <- function(data) {
    if (all(data == "FALSE")) {
      return("FALSE") 
    } else if (all(data %in% c("FALSE", "UP"))) {
      return("UP")
    } else if (all(data %in% c("FALSE", "DOWN"))) {
      return("DOWN")
    } else {
      return(paste("Disagreeing probes:", sum(data == "UP"), "up,", 
                   sum(data == "DOWN"), "down"))
    }
  }
  
  dataframe %>%
    group_by(Day, Concentration, TargetID) %>%
    summarize(
      ENTREZ_GENE_ID = ENTREZ_GENE_ID[1],
      logFC = absMax(logFC),
      #Symmetrical_raw_fold_change = absMax(Symmetrical_raw_fold_change),
      adj.P.Val = min(adj.P.Val),
      Direction = getDir(Direction)
    )
}

#This is a data frame that has one row per GENE (not probe) and the logFC 
# and p val have been summarized as the value for which ever of the probes
#for that gene has the biggest absolute value. If the probes disagree, it says
#how many are up and down

collapsedProbes<-CLcollapseProbes(processedData)

#separate the data with no disagreements from the disagree-ing

collapsedProbesAgree<-filter(collapsedProbes,
                             Direction == "DOWN" | Direction == "UP")
collapsedProbesNOTAgree<-filter(collapsedProbes,
                                Direction != "DOWN" & Direction != "UP")

#plot like it's hot
geneCountAgreePlot<-ggplot(collapsedProbesAgree,aes(x = Day))+
  scale_fill_manual(values = c("UP"="coral", "DOWN"= "cornflowerblue"))+
  geom_bar(aes(fill = Direction),position="dodge")+
  facet_wrap(~Concentration)+
  ggtitle("Count of up and down regulated genes:only agreeing probes")

ggsave("geneCountAgreePlot.png",dpi=600)

geneCountNOTAgreePlot<-ggplot(collapsedProbesNOTAgree,aes(x = Day))+
  geom_bar(aes(fill = Direction),position="dodge")+
  facet_wrap(~Concentration)+
  ggtitle("Count of disaggreeing probes")

ggsave("geneCountNOTAgreePlot.png",dpi=600)
###################COMPARING WITH CYBER T RESULTS ###########################################

#Lamar's cyber T results 
load("LMF CyberT Analysis.Rda")

#Where do limma nonfalse probes overlap with CyberT nonfalse probes?



######################## EXTRACT AND SUBSET CYBERT NON FALSE #################
#function for subsetting nonFalse CyberT data
library(stringr)
library(biomaRt)
selectCyberT<-function(day.dose){
  x<-CyberT[,str_detect(colnames(CyberT),day.dose)]#extract using day and dose
  
  y<-cbind(x[,3:4],CyberT$Probe.ID)#include just pval and direction,add the probe column
  
  names(y)<-c("PVal","Direction","Probe.ID") #change names so all dfs have same
  
  y[y$Direction!=FALSE,]#remove FALSE probes
}
#list of the parts of the column names I want to capture with str_detect

dayDoseList<-list("D1.50v","D1.500","D4.50v","D4.500","D7.50v","D7.500",
                  "D14.50v","D14.500")

#apply the dayDose list with the function
CyberTDfs<-lapply(dayDoseList, FUN=selectCyberT)

#how many non-falses are there?
lapply(CyberTDfs,FUN=nrow)






















############################ FIND OVERLAPS ##########################
findOverlaps<-function(df){
  C<-CyberTDfs[[df]]
  L<-PValLogFCfiltered[[df]]
 list(C[C$Probe.ID %in% L$Probe.ID,],L[L$Probe.ID %in% C$Probe.ID,])
}


#overlaps is a list with 8 elements one for each data group
# each data group has 2 elements, subset of Cyber that overlaps with Limma
# and subset of limma that overlaps with Cyber
overlaps<-lapply(GroupList,FUN = findOverlaps)

#extract the first element of the 8th element
overlaps[[8]][[1]]





######################## LOOKING AT GROUP 8 ONLY ##############
limma8<-overlaps[[8]][[2]]
cyber8<-overlaps[[8]][[1]]


overlapLimma<-limma8%>%
  select(Probe.ID,adj.P.Val)%>%
  mutate(Analysis = rep("Limma", times = nrow(limma8)))%>%
  arrange(adj.P.Val)%>%
  mutate(rankingL = seq(1:nrow(limma8)))%>%
  arrange(Probe.ID)

colnames(overlapLimma)<-c("Probe","Pval","Analysis","ranking")

overlapCyber<-cyber8%>%
  mutate(Analysis = rep("CyberT", times = nrow(cyber8)))%>%
  arrange(PVal)%>%
  mutate(rankingC = seq(1:nrow(cyber8)))%>%
  arrange(Probe.ID)%>%
  select(Probe.ID,PVal,Analysis,rankingC)
  

colnames(overlapCyber)<-c("Probe","Pval","Analysis","ranking")

together<-rbind(overlapCyber,overlapLimma)

together$Probe<-as.character(together$Probe)

ggplot(together, aes(Probe, Pval))+
  geom_point(aes(color = Analysis, size=ranking))





#What about the difference between the ranks?

diff<-data.frame("rankDiff"=overlapCyber$ranking - overlapLimma$ranking,
                 "Probe" = overlapLimma$Probe)

diff$Probe<-as.character(diff$Probe)


ggplot(diff,aes(Probe,rankDiff))+
  geom_point(aes(), size=4, alpha=0.2)
