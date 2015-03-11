
require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
require(pander)

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes")

dataSummary<-read.csv("ex_vivo_HVE_Tenofovir_Summary_results.csv")


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

ggplot(meltedupDownCount,aes(x=day.dose,y=Count))+
  geom_point(aes(color=Direction),size=4)+
  geom_line(aes(group=Direction))

##I am going to try to use all the probes on the array 
#as a background reference for annotation and functional
#analysis stuff. So I will try to extract that info from
# the raw data file that LMF sent

#can't read in data using lumiR, says there is nothing for the
#control data slot...

#this tells me that HumanHT12v4_121001 was the assay used
#require(lumi)
#x<-lumiR("FinalReport_exvivo_TNF_HVE.txt")
#readLines("FinalReport_exvivo_TNF_HVE.txt", n=10)
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


#using the function to make data in longform...
source("DEG-to-long-form.R")

#a function to extract and write data from specific days and concentrations
#the table is for innateDB.
#the csv is for david where you just put in the IDs
#For innateDB upload: first col is entrez, second is log2foldchange, 3rd is adj pval)
extract<-function(day,concentration,direction){
x<-longForm[longForm$Day == day &
              
              longForm$Concentration == concentration &
              
              longForm$Direction == direction,
            c(3,6,8)]

write.table(x,file=paste(day,concentration,direction,"txt",sep="."),
         sep="\t",row.names=FALSE,col.names=FALSE)

write.csv(x[,1],file=paste(day,concentration,direction,"csv",sep="."),
          row.names=FALSE)
}

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes")
#all the UPs from concentration=50
extract("1","50","UP")
extract("4","50","UP")
extract("7","50","UP")
extract("14","50","UP")


#get DAVID data (saved as tab delim txt)
#left the defaults but chose "Fisher's exact" under options
# and chose to include a column for Foldchange

#here is function for getting the DAVID data out of the folder
#I also added a column to each df saying which day it represented
getDAVID<-function(day, concentration, direction){
 y<-read.csv(file=paste("DAVID",day, concentration,direction,"csv", sep="."),
               header=TRUE)
}

#extract the david data from the folder where I saved it
DAVID.1.50.UP<-getDAVID("1","50","UP")
DAVID.1.50.UP<-DAVID.1.50.UP%>%
  mutate(Day= rep("1", times=nrow(DAVID.1.50.UP)))


DAVID.4.50.UP<-getDAVID("4","50","UP")
DAVID.4.50.UP<-DAVID.4.50.UP%>%
  mutate(Day= rep("4", times=nrow(DAVID.4.50.UP)))


DAVID.7.50.UP<-getDAVID("7","50","UP")
DAVID.7.50.UP<-DAVID.7.50.UP%>%
  mutate(Day= rep("7", times=nrow(DAVID.7.50.UP)))


DAVID.14.50.UP<-getDAVID("14","50","UP")
DAVID.14.50.UP<-DAVID.14.50.UP%>%
  mutate(Day = rep("14", times=nrow(DAVID.14.50.UP)))

#combine all the data frames
allDAVID.50.UP<-rbind(DAVID.1.50.UP,DAVID.4.50.UP,DAVID.7.50.UP,DAVID.14.50.UP)

#make Day a factor
allDAVID.50.UP$Day<-factor(allDAVID.50.UP$Day, levels=c("1","4","7","14"))
                           

##DF showing overlapping terms for dose = 50 days 1,4,7,14
overlap.DAVID.50.UP<-allDAVID.50.UP%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


#if you do summarize(n())where filter is  it shows you the term and 
#how many occurances there were in overlapTerms.50.UP. We only want the terms where
#there there 4 occurances (1 per day we looked at)
#a short version for easy viewing

pander(overlap.DAVID.50.UP %>%
  select(Day,Term,PValue,Benjamini),file = "overlap.50.UP")

write.csv(overlap.DAVID.50.UP %>%
            select(Day,Category,Term,PValue,Benjamini),
          file = "overlap.DAVID.50.UP.csv")


###plot plot plot
#note that I changed the x axis labels to a shortened version
# of the term for readability

#Pvalues
ggplot(data=overlap.50.UP, aes())+
  geom_point(aes(x = Term , y =PValue,color=Day),
             position=position_jitter(w=0.15),size=4)+
  geom_hline(y=0.05)+
  scale_x_discrete(labels=c("Acetylation","Cytoskeleton",
                            "Nuclear body","Nuclear speck"))+
  ggtitle("GSEA P-values for overlapping GO terms \n\ in Up-regulated genes,dose=50")
  

#Fold enrichment
ggplot(data=overlap.50.UP, aes())+
  geom_point(aes(x = Term , y =Fold.Enrichment,color=Day),
             position=position_jitter(w=0.15),size=4)+
  scale_x_discrete(labels=c("Acetylation","Cytoskeleton",
                            "Nuclear body","Nuclear speck"))+
  ggtitle("Fold Enrichment for overlapping GO terms \n\ in Up-regulated genes,dose=50")

#Count of genes relating to the term

ggplot(data=overlap.50.UP, aes())+
  geom_point(aes(x = Term , y =Count,color=Day),
             position=position_jitter(w=0.15),size=4)+
  scale_x_discrete(labels=c("Acetylation","Cytoskeleton",
                            "Nuclear body","Nuclear speck"))+
  ggtitle("Number of up-regulated genes associated with overlapping terms \n\
dose=50")

#################Overlapping terms not incl day 1###########\

allDAVID.not1.50.UP<-allDAVID.50.UP%>%
  filter(Day !=1)


overlap.DAVID.not1.50.UP<-allDAVID.not1.50.UP%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 3 days
  arrange(Day,Term,PValue)%>%
  ungroup()


pander(overlap.DAVID.not1.50.UP %>%
         select(Day,Term,PValue,Benjamini))

write.csv(overlap.DAVID.not1.50.UP %>%
            select(Day,Category,Term,PValue,Benjamini),
          file = "overlap.DAVID.not1.50.UP.csv")

#pvalues for those overlapping terms, more low ones for d7
ggplot(data=overlap.not1.50.UP, aes())+
  geom_point(aes(x = Day , y = PValue),size=4,alpha=0.1)+
  geom_hline(y=0.05)+
  ggtitle("Uncorrected Pvalues of overlapping upreg terms in day 4,7,& 14
          \n\ dose = 50")

###############Dose = 50, DOWN ##########################
#extracting data from the long form summary
extract("1","50","DOWN")
extract("4","50","DOWN")
extract("7","50","DOWN")
extract("14","50","DOWN")

#Now read in the DAVID data
DAVID.1.50.DOWN<-getDAVID("1","50","DOWN")
DAVID.1.50.DOWN<-DAVID.1.50.DOWN%>%
  mutate(Day= rep("1", times=nrow(DAVID.1.50.DOWN)))

DAVID.4.50.DOWN<-getDAVID("4","50","DOWN")
DAVID.4.50.DOWN<-DAVID.4.50.DOWN%>%
  mutate(Day= rep("4", times=nrow(DAVID.4.50.DOWN)))

DAVID.7.50.DOWN<-getDAVID("7","50","DOWN")
DAVID.7.50.DOWN<-DAVID.7.50.DOWN%>%
  mutate(Day= rep("7", times=nrow(DAVID.7.50.DOWN)))

DAVID.14.50.DOWN<-getDAVID("14","50","DOWN")
DAVID.14.50.DOWN<-DAVID.14.50.DOWN%>%
  mutate(Day= rep("14", times=nrow(DAVID.14.50.DOWN)))

allDAVID.50.DOWN<-rbind(DAVID.1.50.DOWN,DAVID.14.50.DOWN,
                        DAVID.4.50.DOWN,DAVID.7.50.DOWN)

allDAVID.50.DOWN$Day<-factor(allDAVID.50.DOWN$Day,
                             levels=c("1","4","7","14"))

overlap.DAVID.50.DOWN<-allDAVID.50.DOWN%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()

pander(select(overlap.DAVID.50.DOWN,Day,Term,PValue,Benjamini))

write.csv(overlap.DAVID.50.DOWN %>%
            select(Day,Category,Term,PValue,Benjamini),
          file = "overlap.DAVID.50.DOWN.csv")
############### overlap 50 DOWN not incl day 1##############

allDAVID.not1.50.DOWN<-allDAVID.50.DOWN%>%
  filter(Day !=1)


overlap.DAVID.not1.50.DOWN<-allDAVID.not1.50.DOWN%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 3 days
  arrange(Day,Term,PValue)%>%
  ungroup()


pander(overlap.DAVID.not1.50.UP %>%
         select(Day,Term,PValue,Benjamini))


write.csv(overlap.DAVID.not1.50.DOWN %>%
            select(Day,Category,Term,PValue,Benjamini),
          file = "overlap.DAVID.not1.50.DOWN.csv")

########################Innate DB data#######################

#The data from the browser didn't show any significant pvalues
#they all appear to be 1...I saved  xls, saved as csv and then
#read in with read.csv. Txt files gave problems with some rows
#having missing data. and the "number of items is not a multiple of the number of columns"
#warning

getInnate<-function(day, concentration, direction){
read.csv(file = paste("Innate",day, concentration,direction,"csv", sep="."),
           header=TRUE)
}

#there is a better way to do this with an apply function
#but I haven't figured it out yet.
Innate.1.50.UP<-getInnate("1","50","UP")

Innate.4.50.UP<-getInnate("4","50","UP")

Innate.7.50.UP<-getInnate("7","50","UP")

Innate.14.50.UP<-getInnate("14","50","UP")

Innate.14.50.UP<-getInnate("14","50","UP")

#a function to find instances where pvalue !=1

not1<-function(df){
df[df$Pathway.p.value..corrected.!= 1,]
}
 
# a list of the dfs
dfList<-list(Innate.1.50.UP,Innate.4.50.UP,
               Innate.7.50.UP,Innate.14.50.UP)

#apply not1 over the list
lapply(dfList,FUN=not1)

#looks like there aren't any rows where the pvalue !=1 except
#in the day 7 df where there is a row with NAs

##now the down innate data for dose = 50

Innate.1.50.DOWN<-getInnate("1","50","DOWN")
Innate.4.50.DOWN<-getInnate("4","50","DOWN")
Innate.7.50.DOWN<-getInnate("7","50","DOWN")
Innate.14.50.DOWN<-getInnate("14","50","DOWN")

dfList2<-list(Innate.1.50.DOWN,Innate.4.50.DOWN,
             Innate.7.50.DOWN,Innate.14.50.DOWN)

lapply(dfList2, FUN=not1)
sum(is.na(Innate.7.50.DOWN))

sum(is.na(Innate.7.50.UP))
identical(Innate.7.50.DOWN,Innate.7.50.UP)
#looks like 5 rows with NA in both UP and DOWN for 7.50
#but they aren't the same so I didn't read in the wrong thing


############ biomaRt#############################
#Goal: get GO ids and name from bioMaRt from list of entrez ids
UP.50.1<-read.table("1.50.UP.txt", sep="\t", header = FALSE)
head(UP.50.1)
#make a vector of the ids to annotate
entrez<-as.vector(UP.50.1$V1)

#set the mart you want to use (see listMarts())
#and choose the dataset to use see listDatasets()

ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#the attribute for go term name is "name_1066"
goids<-getBM(attributes=c("entrezgene","go_id","name_1006"),
            filters = "entrezgene",values=entrez,
            mart = ensembl)

