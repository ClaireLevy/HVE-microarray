#Working with raw illumina data
#started 26March15


require(dplyr)
library(lumi)
library(limma)
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData")
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

#select all the such that the columns arranged by dose and then day
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
#To do that, I pass this logical vector to the lumiBatch to
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
#Just looking at day14 dose 500 (and controls)
load(dataMatrixSDfilter)
#select columns with day 14 data 
dose500day14Data<-dataMatrixSDfilter[,c("HVE_A4", "HVE_A8","HVE_A12",
                                   "HVE_C4", "HVE_C8","HVE_C12")]

#making a design matrix based on limma vignette pg 42


#set up the targets frame, making sure the sample order
#is the same as in the data (3 controls, 3 treatments)
targets<-data.frame("CellLine"=c("HVE1","HVE2","HVE3",
                                 "HVE1","HVE2","HVE3"),
                    "Treatment"=c("control","control","control",
                                  "drug","drug","drug"))

#using the targets frame, make a design matrix that accounts for
#pairing within Cell lines (i.e. HVE1 dose=0 and HVE1 dose=500 
#go together). see limma guide pg 42

CellLine<-factor(targets$CellLine, levels=c("HVE1","HVE2","HVE3"))
Treatment<-factor(targets$Treatment, levels=c("control","drug"))

#make the design matrix
design<-model.matrix(~0 + CellLine + Treatment)

#Do I need a contrasts matrix?

#fit the model using the design matrix and the data
fit<-lmFit(dose500day14Data,design)

#"Given a microarray linear model fit, compute moderated t-statistics,
#moderated F-statistic, and log-odds of differential expression by 
#empirical Bayes moderation of the standard errors towards
#a common value."

fit<-eBayes(fit)

topTable(fit,coef="Treatmentdrug",adjust="BH")

#Now need to filter for logFC using decideTests??



