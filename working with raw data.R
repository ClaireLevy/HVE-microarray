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


summary(RAW.lumi,"QC")
#some QC plots

preNormDensity<-density(RAW.lumi)
#looks like HVE_A1 is lower than most others

preNormRevCDF<-plotCDF(RAW.lumi,reverse=TRUE)
#CDF is the probability of a getting an intensity <= to x
#A PDF answers the question: "How common are samples at exactly this value?" 
#A CDF answers the question "How common are samples that are less than or equal to this value?"
#The CDF is the integral of the PDF.


preNormBox<-boxplot(RAW.lumi)
#this doesn't really mean anything to me except that some cluster
#together and others don't.
plotSampleRelation(RAW.lumi, method="mds")

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



#check the plots after normalization and transformation
#all more uniform
NormDensity<-density(RAWlumi.N.Q)
NormRevCDF<-plotCDF(RAWlumi.N.Q, reverse=TRUE)
NormBox<-boxplot(RAWlumi.N.Q)

plotSampleRelation(RAWlumi.N.Q, method="mds")

head(detection(RAWlumi.N.Q))

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
  
#test this code by making up a matrix where some rows meet the condition
#that 3 consecutive columns are all<0.05 and some dont.
test<-read.csv("testMatrix.csv")
test
#The resulting vector should be T, T, T ,F,T,F,F
print(testLogical<-rowSums(test[,1:3]<0.05)==3|
  rowSums(test[,4:6]<0.05)==3|
  rowSums(test[,7:9]<0.05)==3)
#looks ok



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
dataMatrix<-exprs(expressed)
hist(dataMatrix, breaks=50)

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
load("dataMatrixSDfilter.Rda")
library(limma)
library(genefilter)


#sample key has sample information incl a "Group" column
#which assigns a number to a day + treatment + paired control (8)
#ex) day 1+dose 50+day1 no drug = group1
#day 1 is groups 1&2 (dose 50 then 500), day 2 is groups 3&4 etc


SampleKey<-read.csv("SampleKey.csv",stringsAsFactors=FALSE)

#Function to select a particular group from the data set

SelectDataGroup<-function(Group){
  SampleKey[SampleKey$Group == Group,]
  dataMatrixSDfilter[,x$SampleName]
}

# I want a list consisting of data frames for each group
#Here is a list of the groups
GroupList<-list(1,2,3,4,5,6,7,8)
#pass the list to the SelectDataGroup function
GroupDataframes<-lapply(GroupList,FUN=SelectDataGroup)


#This file has the code for making the design matrix and targets df
source("DesignMatrixAndTargets.R")

#fit the model using the design matrix and the data selected above
fit<-lapply(GroupDataframes,FUN=lmFit,design = design)


#"Given a microarray linear model fit, compute moderated t-statistics,
#moderated F-statistic, and log-odds of differential expression by 
#empirical Bayes moderation of the standard errors towards
#a common value."

fit<-lapply(fit, FUN=eBayes)

#generate a top table from the treatment data with BH adj pvalues
#give all entries (so as many rows as are in fit)

TT<-lapply(fit,FUN=topTable,coef="Treatmentdrug",
           adjust="BH", number = 10380)

#filter for logFC
LogFCfiltered<-lapply(TT,subset,logFC>=0.5 | logFC<=-0.5 )
#filter for adj p value
PValLogFCfiltered<-lapply(filtered,subset,adj.P.Val <=0.05)

##############################################################
#Lamar's cyber T results (after doing ORA)
CyberT<-read.csv("J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\differentiallyExpressedGenes\\ex_vivo_HVE_Tenofovir_Summary_results.csv")

save(CyberT, file = "LMF CyberT Analysis.Rda")

filteredTT500day14$Probe<-as.integer(filteredTT500day14$Probe)
CyberT$Probe.ID<-as.integer(CyberT$Probe.ID)


#any overlaps?
sum(CyberT$Probe.ID %in%)
#62

#First remove all but the summary columns and probes that overlap


######################### You need to do long vector %in% shorter vector or else you get
######################### a zillion rows for some reason
# sum(filteredTT500day14$Probe %in% CyberT$Probe.ID) = 62 BUT
### This: CyberT[filteredTT500day14$Probe %in% CyberT$Probe.ID,c(1,40)]
##### Give you 6838 observations, not 62

#subset for just the overlapping probes and select probe.id and day 14 dose 500
CyberTshort<-CyberT[CyberT$Probe.ID%in% filteredTT500day14$Probe,
                    c(1,40)]

#so 62 of my 71 probes (not ORA analyzed) were present in Cyber T analysis

#remove "FALSE" entries (don't meet logFC and pval cutoffs)
NoFalseCyberTshort<-filter(CyberTshort,D14.500vNT.DEG.1!="FALSE")


#How many probes that passed cutoffs in my analysis (no ORA yet)
#were also present AND passed cutoffs in CyberT?
nrow(NoFalseCyberTshort)

#so, comparing what I found for day14 dose 500 to cyberT findings:

overlap_summary<-data.frame("NProbes"=nrow(topTable500day14),
                            "meet_cutoffs"=nrow(filteredTT500day14),
                            "present_in_Cyber_T"=nrow(CyberTshort),
                            "meet_cutoffs_in_Cyber_T"=nrow(NoFalseCyberTshort)
)
###############################try the same thing with day14 dose= 50

#select columns with day 14 dose = 50 data 
dose50day14Data<-dataMatrixSDfilter[,c("HVE_A4", "HVE_A8","HVE_A12",
                                        "HVE_B4", "HVE_B8","HVE_B12")]


#fit the model using the design matrix and the data
fit<-lmFit(dose50day14Data,design)
fit<-eBayes(fit)

#generate a top table from the treatment data with BH adj pvalues
#give all entries (so as many rows as are in fit)
topTable50day14<-topTable(fit,coef="Treatmentdrug",adjust="BH", number = nrow(fit))


#make the probe ids an actual column, not just the rownames
topTable50day14$Probe<-rownames(topTable50day14)
#remove the rownames "column"
rownames(topTable50day14)<-NULL


#filter for |logFC|>=0.5 and adj p <=0.05
filteredTT50day14<-topTable50day14%>%
  filter(logFC>=0.5 | logFC<=-0.5)%>%
  filter(adj.P.Val <= 0.05)

CyberTshort2<-CyberT[CyberT$Probe.ID%in% filteredTT50day14$Probe,
                    c(1,36)]


#remove "FALSE" entries (don't meet logFC and pval cutoffs)
NoFalseCyberTshort2<-filter(CyberTshort2,D14.50vNT.DEG!="FALSE")

