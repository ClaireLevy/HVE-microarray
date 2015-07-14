####################### EXTRACT DATA FOR CYBER T INPUT ###############
setwd("J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\")
library(dplyr)
# need to submit data to cyber T webform like so:
# Label, Cond1Repl1, Cond1Repl2, Cond2Repl1, Cond2Repl2
#dataMatrixSDfilter is already in this order
#(HVE1cntl, HVE2cntl, HVE3cntl, HVE1treat, HVE2treat, HVE3treat)
#but there isn't a labled for the probe column so I need to make
#it a column instead of rownames
load("dataMatrixSDfilter.Rda")


load("SampleKey.Rda")

#Here is a vector of the groups
GroupList<-list(1,2,3,4,5,6,7,8)

#Function to select a particular group from the limma data set
getForCyber<-function(Group){
   x<-SampleKey[SampleKey$Group == Group,]
   m<-as.data.frame(dataMatrixSDfilter[,x$SampleName])
   m$Probe<-rownames(m)
   rownames(m)<-NULL
   m<-m[ ,c(7,1:6)]
   write.csv(m, paste("Group",Group,"csv",sep="."),row.names=FALSE)
   m
}

#pass the list to the getForCyber function

p<-lapply(GroupList,FUN=getForCyber)
#check that my function works ok
t<-getForCyber(1)#using function manually (no apply)
t2<-p[[1]]#using apply


#extract group 1 data without using a function
t3<-as.data.frame(dataMatrixSDfilter[,1:6])
t3$Probe<-rownames(t3)
rownames(t3)<-NULL
t3<-t3[,c(7,1:6)]

identical(t,t2) #true
identical(t2,t3)#true

# structure of example input file
head(read.csv("Group.1.csv"))
# Probe   HVE_A1    HVE_A5   HVE_A9   HVE_B1    HVE_B5   HVE_B9
# 4760377 8.311434  7.931828 8.239210 8.233480  8.019671 8.329865
# 4050154 9.528693 11.097650 7.634868 9.256283 10.654425 7.602925
# 3610072 8.186686  9.080132 8.397511 8.363603  8.657692 8.300633
#   60689 8.008220  8.345514 8.468531 8.064208  8.302465 8.536110
# 1580504 8.198067  8.302021 7.859298 8.034230  8.290039 7.962794
# 7650020 8.217213  8.163575 8.470606 8.349379  8.139904 8.314603

#in cyberT http://cybert.ics.uci.edu/pair/
#normalization = NONE
#sliding window size = 101
#Bayesian confidence value = 10
#no post processing
####################### READ IN CYBER T RESULTS ####################
#Cyber T output columns:
# Lab_1 : Label column input by user
# R_1 : Ratio column #1 input by user
# R_2 : Ratio column #2 input by user
# nR : Number of ratios

#There is a column here called EstExpr but CyberT doesn't say what it means

# meanR : Mean of the log transformed ratios
# stdR : Standard devation of the ratios
# rasdR : the background standard deviation for ratios (if Bayesian analysis is performed)
# bayesSD : The Bayesian or regularized standard deviation of the ratios (if Bayesian analysis is performed)
# ttest :: The t-test statistic calculated from ratios data using either stdR or bayesSD
# DF : The degrees of freedom for the t-test statistic (nR-1) (if Bayesian analysis is NOT performed)
# bayesDF : The degrees of freedom for the t-test statistic plus that associated with the Bayesian estimate (if Bayesian analysis is performed)
# pVal : The p-value associated with the t-test on ratios (column ttest) with DF or bayesDF degrees of freedom
# Bonferroni : Bonferroni corrected q-values (if multiple hypothesis testing correction is performed)
# BH : Benjamini & Hochberg corrected q-values (if multiple hypothesis testing correction is performed

data <- lapply(1:8, function(n) {
   data <- read.table(paste("Group", n, "cyber.txt", sep = "."), 
      sep = "\t", header = TRUE)
   data$Group <- n
   data
})

data <- do.call("rbind", data)

groupKey <- SampleKey %>% 
   group_by(Group) %>% 
   summarize(Day = Day[1], 
      Concentration = paste(Concentration[1], Concentration[4])) 
data <- merge(data, groupKey)

# Question: is the logFC calculated by Cyber-T by us the same as the logFC 
# calculated by limma? To answer, get limma data and then compare with 
# Cyber-T data.
limmaAnalysis <- lapply(TT, FUN = function(x) {
   x$Probe <- rownames(x) 
   x
})
limmaAnalysis <- do.call("rbind", limmaAnalysis)
limmaAnalysis$Group <- rep(1:8, each = 10380)
limmaAnalysis <- select(limmaAnalysis, logFC, Group, Probe)
cyberAnalysis <- mutate(data, Probe = Lab_0) %>% select(meanR, Group, Probe)

merged <- merge(limmaAnalysis, cyberAnalysis, by = c("Probe", "Group"))
sum(abs(merged$logFC - merged$meanR) < 0.00001)
# [1] 83040
# i.e. yes, for every probe, limma and cyber-t calculate the same log FC

# What about the log fold change from the analysis from the eLife paper?
source("DEG-to-long-form.R")
previous <- mutate(longForm, Probe = Probe.ID, 
   Concentration = paste(0, Concentration)) %>% 
   merge(groupKey) %>% 
   select(Log2_fold_change, Group, Probe)
supermerged <- merge(merged, previous)
sum(abs(supermerged$logFC - supermerged$Log2_fold_change) < 0.00001)
# [1] 2
# The log fold changes are different between what we've done and the analysis
# for the eLife paper, suggesting that we're somehow transforming (or not 
# transforming) the data in a way that's different or that the data file we're
# working from is somehow different. 

# List of things that might explain this: 
# - Different raw data / format
# - Something incorrect in experiment design
# - Missing or extra data transformation

# Summary of Cyber-T results
# no logFC filter
data %>% 
   group_by(Group, Day, Concentration) %>% 
   summarize(sigPVal = sum(pVal < 0.05), sigBH = sum(BH < 0.05), 
      sigBon = sum(Bonferroni < 0.05))
#   Group Day Concentration sigPVal sigBH sigBon
# 1     1   1          0 50     558     0      0
# 2     2   1         0 500     166     0      0
# 3     3   4          0 50    1248     5      2
# 4     4   4         0 500    2218   140      2
# 5     5   7          0 50    2892   378      3
# 6     6   7         0 500    1471    28      1
# 7     7  14          0 50    1743    59      3
# 8     8  14         0 500    2243   514     17

# with logFC filter
data %>% 
   group_by(Group, Day, Concentration) %>% 
   filter(abs(meanR) >= 0.5) %>% 
   summarize(sigPVal = sum(pVal < 0.05), sigBH = sum(BH < 0.05), 
      sigBon = sum(Bonferroni < 0.05))
# Group Day Concentration sigPVal sigBH sigBon
#     1   1          0 50     262     0      0
#     2   1         0 500      11     0      0
#     3   4          0 50     195     5      2
#     4   4         0 500     234   106      2
#     5   7          0 50     815   348      3
#     6   7         0 500     182    27      1
#     7  14          0 50      68    34      3
#     8  14         0 500      71    71     17

# limma 
# group  sig
# 1      0  
# 2      0  
# 3      0  
# 4      0  
# 5      0  
# 6      0 
# 7      15 
# 8      71