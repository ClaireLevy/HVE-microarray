####################### EXTRACT DATA FOR CYBER T INPUT ###############
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

#in cyberT http://cybert.ics.uci.edu/pair/
#normalization = NONE
#sliding window side = 101
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

   
x <- lapply(1:8, function(n) {
   data <- read.table(paste("Group", n, "cyber.txt", sep = "."), 
      sep = "\t", header = TRUE)
   data$Group <- n
   data
})

data <- do.call("rbind", x)

data %>% 
   group_by(Group) %>% 
   filter(abs(meanR) >= 0.5) %>% 
   summarize(sigPVal = sum(pVal < 0.05), sigBH = sum(BH < 0.05), 
      sigBon = sum(Bonferroni < 0.05))
