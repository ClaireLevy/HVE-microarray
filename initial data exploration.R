
require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
dataSummary<-read.csv("ex_vivo_HVE_Tenofovir_Summary_results.csv")
head(dataSummary)

#spread sheet is laid out with the 4 columns: log2 fold change,
#raw fold change, FDR, Differentially expressed gene (DEG)
#These 4 columns are repeated for the 4 time points and for the
#two doses of tenofovir (500 and 50)

#timepoints are D1, D4, D7, D14



#I am going to isolate just the DEG summary data from the end
#of the spreadsheet (columns 36-43) and the annotations in the first
#three columns and then
#extract the probe, target and Entrez IDs for the
#DEG for each time point in each direction (up or down)


summaryOnly<-dataSummary[,c(1:3,36:43)]
head(summaryOnly)
names(summaryOnly)[4:11]<-c("d14.50","d7.50","d4.50","d1.50",
                            "d14.500","d7.500","d4.500","d1.500")



##melt and show HOW MANY are up, down and false for each sample
melted<-melt(summaryOnly,
             id.vars=c("Probe.ID","TargetID","ENTREZ_GENE_ID"),
             variable.name="day.dose",value.name="UpDownFalse")

#change order of samples for future plots
melted$day.dose<-factor(melted$day.dose,
                        levels=c("d1.50","d4.50","d7.50",
                                 "d14.50","d1.500","d4.500","d7.500",
                                 "d14.500"))

upDownCount<-melted%>%
  group_by(day.dose)%>%
  summarise(countUP=sum(UpDownFalse=="UP"),
            countDOWN=sum(UpDownFalse=="DOWN"),
            countFALSE=sum(UpDownFalse=="FALSE"))

require(pander)
pander(upDownCount)

#now melt that df to prepare it for easy plotting
meltedupDownCount<-melt(upDownCount, id.vars="day.dose",
                        variable.name="Direction",value.name="Count")


#how many genes are there total?
nGenes<-length(summaryOnly$Probe.ID)


##plot of up down and false counts
ggplot(meltedupDownCount,aes(x=day.dose,y=Count))+
  geom_bar(aes(fill=Direction), position="dodge",stat="identity")+ylab("Number of genes")+
  geom_hline(y=7832)+
  ggtitle("Number of upregulated, downregulated and non-significant genes\n\ per time point and dose\n\ 
          (line represents total genes analyzed for diff exp)")

  
##I am going to try to use all the probes on the array 
#as a background reference for annotation and functional
#analysis stuff. So I will try to extract that info from
# the raw data file that LMF sent

#can't read in data using lumiR, says there is nothing for the
#control data slot...

#this tells me that HumanHT12v4_121001 was the assay used
require(lumi)
x<-lumiR("FinalReport_exvivo_TNF_HVE.txt")
readLines("FinalReport_exvivo_TNF_HVE.txt", n=10)
#But DAVID doesn't have that as one of their defaults
#I will try to extract them.

# data<-readBeadSummaryData("FinalReport_exvivo_TNF_HVE.txt")
# str(data)
# #looks like featureData slot has columns for ProbeID, TargetID,
# #and PROBE_ID. The last one is illumina (starts with ILMN_)so
# #I am guessing that ProbeID is the same as Probe.ID in the summarydata
# 
# 
# 
# 
# 
# 
# 
# #If ProbeID = Probe.ID (from summary data), there should be some matches if I do
# # %in%...
# str(melted)#they are integers here but factor in the raw data
# #so I'll change them
# 
# melted$Probe.ID<-factor(melted$Probe.ID)
# 
# #now check for matches
# x<-sum(melted$Probe.ID %in% y)
# length(melted$Probe.ID)
# #ok looks like all Probe.IDs are in ProbeID
# 
# 
# z<-featureNames(data)#this apparently gives the feature names
# #for the chip, which I guess is what I want
# y<-featureData(data)$ProbeID
# #what is the difference here? 
# identical(z,y) # I guess they are the same thing
# 
# str(z)
# head(z)
# tail(z)
# 
# str(y)
# head(y)
# tail(y)
# #I don't know why but there are sample names at the end of these
# #and only 47,323 probe Ids.
# #checked 27Feb14, this is right, there should be 47323.
# 
# 
# ## DAVID couldn't figure out the IDs so I'll
# #convert entrez>symbol and use symbol from the raw data
# #nevermind, DAVID won't let you use gene symbol for a background.
# 
# 
# #so I'll try the illumina ids
# 
# str(featureData(data))
# str(data)
# allGenesIlmn<-(featureData(data)$PROBE_ID)
# #remove blanks and NAs
# allGenesIlmn<-allGenesIlmn[allGenesIlmn!=""]
# allGenesIlmn<-na.omit(allGenesIlmn)
# #write
# 
# write.csv(allGenesIlmn[1:47323],"allGenesIlmn.csv")

#DAVID didn't understand all of them so I chose the recommended
#"map the IDs that DAVID could convert" option

d4.50<-melted%>%
  filter(day.dose=="d4.50",UpDownFalse!=FALSE)%>%
  select(ENTREZ_GENE_ID)
d4.50<-na.omit(d4.50)
write.csv(d4.50,"d4.50regulated.csv")

d4.50UP<-melted%>%
  filter(day.dose=="d4.50",UpDownFalse=="UP")%>%
  select(ENTREZ_GENE_ID)
d4.50UP<-na.omit(d4.50UP)

write.csv(d4.50UP,"d4.50upReg.csv")

#reading in data incl log2FC and adj p values for innate db exploration
allData<-read.csv("CL copy of LMF summary results.csv")

#just for d1.50, regulated genes only, both up and down
d1.50<-select(allData, ENTREZ_GENE_ID, d1.50adjP,d1.50FC,d1.50DEG)
d1.50<-filter(d1.50,d1.50DEG!="FALSE")
write.csv(d1.50,"d1.50.csv")
