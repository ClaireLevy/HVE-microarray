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

#select all the data such that the columns arranged by dose and then day
allDays.N.Q<-RAWlumi.N.Q[,c("HVE_A1","HVE_A5","HVE_A9",
                     "HVE_B1","HVE_B5","HVE_B9",
                     "HVE_C1","HVE_C5","HVE_C9",
                     "HVE_A2","HVE_A6","HVE_A10",
                     "HVE_B2","HVE_B6","HVE_B10",
                     "HVE_C2","HVE_C6","HVE_C10",
                     "HVE_A3","HVE_A7","HVE_A11",
                     "HVE_B3","HVE_B7","HVE_B11",
                     "HVE_C3","HVE_C7","HVE_C11",
                     "HVE_A4","HVE_A8","HVE_A12",
                     "HVE_B4","HVE_B8","HVE_B12",
                     "HVE_C4","HVE_C8","HVE_C12")]



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
#for AT LEAST one of the 9 groups (3 doses x 3 donors) across the
#4 days so I use | to say it can be true for this group OR this group...

# The resulting logical vector is all the probes and the value TRUE or
# FALSE for each of them telling me if at least 1/9 groups fits the
#condition that rowSums==3.


expressedLogical<-rowSums(detection(allDays.N.Q[,1:3]) <0.05)==3|
  rowSums(detection(allDays.N.Q[,4:6]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,7:9]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,10:12]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,13:15]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,16:18]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,19:21]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,22:24]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,25:27]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,28:30]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,31:33]) <0.05) ==3|
  rowSums(detection(allDays.N.Q[,34:36]) <0.05) ==3
  


#Now I need to subset the lumibatch to just contain those TRUE probes.
#To do that, I pass this logical vector to the lumiBatch 
# which subsets the lumiBatch to contain only the probes that are TRUE.
expressed<-allDays.N.Q[expressedLogical,]

dims(RAWlumi.N.Q)
dims(expressed)
#19386 probes left
dims(RAWlumi.N.Q)-dims(expressed)
#removed 27937 probes


############################ FILTER SD  ##################


library(limma)
library(genefilter)


#calculate the row standard deviations
rowSdsMatrix<-rowSds(dataMatrix)

hist(rowSdsMatrix,breaks=50)

#subset dataMatrix to get just the probes where the rowSD is >0.2
dataMatrixSDfilter<-dataMatrix[rowSdsMatrix>=0.2,]


hist(dataMatrixSDfilter, breaks=50)

#compare the filtered and unfiltered lumiBatch

dim(dataMatrix)[1]-dim(dataMatrixSDfilter)[1]
#got rid of 9006 probes
dim(dataMatrixSDfilter)[1]
#10380 left


################ FIT MODEL ################################

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

SampleKey<-read.csv("SampleKey.csv",stringsAsFactors=FALSE)

#Function to select a particular group from the data set
SelectDataGroup<-function(Group){
  x<-SampleKey[SampleKey$Group == Group,]
  dataMatrixSDfilter[,x$SampleName]
}


# I want a list consisting of data frames for each group
#Here is a list of the groups
GroupList<-list(1,2,3,4,5,6,7,8)

#pass the list to the SelectDataGroup function
GroupDataframes<-lapply(GroupList,FUN=SelectDataGroup)




#This file has the code for making the design matrix and targets df
source("DesignMatrixAndTargets.R")

#fit the model using the design matrix and the data frames for each
#group
fit<-lapply(GroupDataframes,FUN=lmFit,design = design)

fit<-lapply(fit, FUN=eBayes)

#generate a top table from the treatment data with BH adj pvalues
#give all entries (so, as many rows as are in fit)
#"inf" means show all entries
TT<-lapply(fit,FUN=topTable,coef="Treatmentdrug",
           adjust="BH", number = "inf")

#filter for logFC
LogFCfiltered<-lapply(TT,subset,logFC>=0.5 | logFC<=-0.5 )
#filter for adj p value
PValLogFCfiltered<-lapply(LogFCfiltered,subset,adj.P.Val <=0.05)



##############################################################
#Lamar's cyber T results (after doing ORA)
load("LMF CyberT Analysis.Rda")

#let's look at group 8 (day 14 dose = 500 and dose =0)

group8TT<-PValLogFCfiltered[[8]]

#make the row names into a column called Probe.limma

group8TT$Probe.limma<-rownames(group8TT)

rownames(group8TT)<-NULL
#put in same format as CyberT
group8TT$Probe.limma<-as.integer(group8TT$Probe)

#how many non-false probes did Lamar find for this day?
#filter for day14 dose 500 non=false and select probes and p vals

NoFalseCyberT<-CyberT%>%
  filter(D14.500vNT.DEG.1!=FALSE)%>%
  select(Probe.ID,FDR.Adj.p.val.D14.500vNT )

length(NoFalseCyberT$Probe.ID)
#4789, so 4727 more than I got.

#subset for just the overlapping probes 
NoFalseCyberTshort<-NoFalseCyberT[NoFalseCyberT$Probe.ID%in% group8TT$Probe.limma,]

group8TTshort<-group8TT[group8TT$Probe.limma %in% NoFalseCyberT$Probe.ID,]
#so 46 of my 71 probes were present in Cyber T 
#analysis for day 14

#so, comparing what I found for day14 dose 500 to cyberT findings:

overlap_summary<-data.frame("Number_of_Probes"=nrow(TT[[8]]),
                            "meet_cutoffs"=nrow(group8TT),
                            "also_found_in CyberT"=nrow(NoFalseCyberTshort)
)

#What about the order of the probes?
#make a df with the data from cyberT and limma that overlaps, arranged
# by p value

overlapLimma<-group8TTshort%>%
  select(Probe.limma,adj.P.Val)%>%
  mutate(Analysis = rep("Limma", times = nrow(group8TTshort)))%>%
  arrange(adj.P.Val)%>%
  mutate(rankingL = seq(1:nrow(group8TTshort)))%>%
  arrange(Probe.limma)

colnames(overlapLimma)<-c("Probe","Pval","Analysis","ranking")
overlapCyberT<-NoFalseCyberTshort%>%
  mutate(Analysis = rep("CyberT", times = nrow(NoFalseCyberTshort)))%>%
  arrange(FDR.Adj.p.val.D14.500vNT)%>%
  mutate(rankingC = seq(1:nrow(NoFalseCyberTshort)))%>%
  arrange(Probe.ID)



colnames(overlapCyberT)<-c("Probe","Pval","Analysis","ranking")




together<-rbind(overlapCyberT,overlapLimma)

together$Probe<-as.character(together$Probe)

ggplot(together, aes(Probe, Pval))+
  geom_point(aes(color = Analysis, size=ranking))

#What about the difference between the ranks?

diff<-data.frame("rankDiff"=overlapCyberT$ranking - overlapLimma$ranking,
                 "Probe" = overlapLimma$Probe)

diff$Probe<-as.character(diff$Probe)


ggplot(diff,aes(Probe,rankDiff))+
  geom_point(aes(), size=4)

#What about the other days/doses?
# subsetting CyberT data

selectCyberT<-function(day.dose){
x<-CyberT[,str_detect(colnames(CyberT),day.dose)]#extract using day and dose
y<-cbind(x[,3:4],CyberT$Probe.ID)#add the probe column
names(y)<-c("PVal","Direction","Probe.ID") #change names so all dfs have same
y[y$Direction!=FALSE,]#remove FALSE probes
}

dayDoseList<-list("D1.50v","D1.500","D4.50v","D4.500","D7.50v","D7.500",
                  "D14.50v","D14.500")
CyberTDfs<-lapply(dayDoseList, FUN=selectCyberT)

nrow(CyberTDfs[[1]])
lapply(CyberTDfs,FUN=nrow)
  
