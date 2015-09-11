#Working with raw illumina data
#started 26March15


require(dplyr)
library(lumi)
library(limma)
library(dplyr)
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


RAW.lumi



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

# from the detection data of the lumiBatch created above,
#select the probes (rows) where detection is <0.05 for ALL of 
#the donors (==3) for dose=0  OR dose = 50 OR dose = 500 on
#any of the days.
#day1 is cols 1:9, day4 is 10:18, etc, 3 donors per dose per day.

#in other words....

#For each row (probe) determine whether the detection p values
#in a group of ALL the donors for a single condition
#(dose + day,ex HVE_A1, A5,A9) are <0.05.
# If they are, rowSums==3 b/c "1"(yes) for all 3 of the columns.
# and that probe is TRUE. If not, FALSE.
#I want to keep probes where the rowSums==3 is TRUE
#for AT LEAST one of the  3 doses x 3 donors across the
#4 days so I use | to say it can be true for this condition OR this condition...

# The resulting logical vector is all the probes and the value TRUE or
# FALSE for each of them telling me if at least 1/9 conditions fits the
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

hist(rowSdsMatrix,breaks=50)

#subset dataMatrix to get just the probes where the rowSD is >0.2
dataMatrixSDfilter<-dataMatrix[rowSdsMatrix>=0.2,]


hist(dataMatrixSDfilter, breaks=50)

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
##GroupDataframes.Rda is a list of dataframes. There is one data frame
# for each sample group. A group is: day + treatment + paired control
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
#give all entries (so, as many rows as are in fit)

TT<-lapply(fit,FUN=topTable,coef="Treatmentdrug",
           adjust="BH", number = Inf)

######################## FILTER FOR LOG FC AND PVAL ###############

#filter for logFC
LogFCfiltered<-lapply(TT,subset,logFC>=0.5 | logFC<=-0.5 )

#filter for adj p value
PValLogFCfiltered<-lapply(LogFCfiltered,subset,adj.P.Val <=0.05)


#change the rowsnames(probes) to an actual column

PValLogFCfiltered<-lapply(GroupList, function(df){
  c<-PValLogFCfiltered[[df]]
  c$Probe.ID<-rownames(c)
  rownames(c)<-NULL
  return(c)
  
})

lapply(PValLogFCfiltered,FUN=nrow)
###################COMPARING WITH CYBER T RESULTS ###########################################

#Lamar's cyber T results 
load("LMF CyberT Analysis.Rda")

#Where do limma nonfasle probes overlap with CyberT nonfalse probes?



######################## EXTRACT AND SUBSET CYBERT NON FALSE #################
#function for subsetting nonFalse CyberT data
library(stringr)
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
#overlaps[[8]][[1]]





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
  geom_point(aes(), size=4)