
require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
require(pander)
source("getAllDavid.R")
source("getAllInnateDB.R")
source("analyzeAndWriteDavid.R")
source('subsetToOverlappingGoTerms.R')

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
############################################################

# To write .txt and .csv files for input into DAVID/InnateDB, use 
# write-files-for-DAVID-and-InnateDB.R

###############DAVID data for dose = 50 UP ######################################
#NOTE: on DAVID site I used the illumina HT 12 v3 background 
# and  in the functional annotation chart,under options
#checked "Fisher exact", then reran with options and saved.

# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

#(df, concentration, filename, direction, day1Sig = FALSE)
analyzeAndWriteDavidWithDay1(allD50, 50, "overlap.DAVID.50.UP.csv", "UP")

###plot plot plot

#Pvalues plot
plot150UP<-ggplot(data=overlapD50up, aes())+
  geom_point(aes(x = Day , y = Benjamini),
             position=position_jitter(w=0.15),
             size=3, 
             alpha = 0.7)+
  geom_hline(y=0.05)+
  ggtitle("P-values for overlapping GO terms \n\ in Up-regulated genes,dose=50")
  

# This will save a 400x400 file at 100 ppi
ggsave("plot150UP.png", width=4, height=4, dpi=100)


#########DAVID data for dose = 50 UP not incl day 1###########\
# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

analyzeAndWriteDavid(allD50, 50, "overlap.DAVID.not1.50.UP.csv", "UP")

###############DAVID data for dose = 50 DOWN ############
# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

analyzeAndWriteDavidWithDay1(allD50, 50, "overlap.DAVID.50.DOWN.csv", "DOWN")

##############DAVID data for dose = 50 DOWN not incl day 1##############
# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

analyzeAndWriteDavid(allD50, 50, "overlap.DAVID.not1.50.DOWN.csv", "DOWN")

############D data for dose = 500 UP #############
# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)

analyzeAndWriteDavidWithDay1(allD500, 500, "overlap.DAVID.500.UP.csv", "UP", TRUE)

###############DAVID data for dose = 500 UP excluding day 1##############
# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)

analyzeAndWriteDavid(allD500, 500, "overlap.DAVID.not1.500.UP.csv", "UP")


###############DAVID data for dose = 500 DOWN ###########
##################DAVID DATA NEEDS TO BE CORRECTED##########

# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)


#####excluding day 1##############

##function to separate out the GO ids in the DAVID data
##so they can be compared with innate data
getGO<-function(df){
  mutate(df,
         Pathway.Id= ifelse(str_detect(df$Term,"GO:")==TRUE,
                            substr(df$Term,1,10),NA))
}



########################Innate DB data#######################

#########Innate data dose = 50###############################
# read in InnateDB data
combinedInnate50<-getAllInnateDB(50)








#which GO ids are in both innate and David
# if the datasets are already filtered for days 
#4,7,14 overlaps?

################DAVID + Innate data (not1) for dose= 50 DOWN

#GO IDs that occur on days 4,7,14 in innate50UP data
overlapInnate50not1UP<-combinedInnate50 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Pathway.Id)%>%
  filter(n()==3)%>%
  ungroup()

# how many GO terms have p < 0.05 on all 3 days, counted by 2 methods
overlapInnate50not1UP %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0)) %>%
   ungroup()  %>% 
   summarize(allSig = sum(allSig))
# says 27 
z <- dcast(mutate(overlapInnate50not1UP, Day = paste0("d", Day)),  Pathway.Id~Day,
   value.var = "Pathway.p.value.corrected")
sum(z$d4 < 0.05 & z$d7 < 0.05 & z$d14 < 0.05)
# says 27 again


#separating go terms out from the pre-existing david data
overlapDnot150up <- subsetToOverlappingGoTerms(getAllDavid(50), "UP", 
   withDay1 = FALSE)
overlapDnot150up<-getGO(overlapDnot150up)

#merge the data sets
ID50UP<-merge(overlapInnate50not1UP,overlapDnot150up,
         by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID50UP<-ID50UP%>%
  select(Day,Pathway.Id,Term,Pathway.p.value.corrected,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value.corrected<0.05)

names(ID50UP)[4:5]= c("InnateDB.p.value","DAVID.p.value")

ID50UP$Day<-factor(ID50UP$Day, levels = c("4","7","14"))


# PLOTS : InnateDB and DAVID p values
ggplot(ID50UP,aes(x = Day, y = value,
                         color = variable))+  
  geom_point(aes(y=InnateDB.p.value,                 
                 col="InnateDB.p.value"), 
             position=position_jitter(w=0.15),size=3,alpha=0.5)+
  geom_point(aes(y=DAVID.p.value,
                 col="DAVID.p.value"),
             position=position_jitter(w=0.15),size=3,alpha=0.5)+
  labs(y="Benjamini")+
  theme(legend.title=element_blank())+
  ggtitle("Overlapping upreg GO terms from InnateDB and DAVID ORA \n\
          dose = 50")



#InnateDB VS DAVID pvalues, limited 0-0.05
ggplot(ID50UP, aes(x = DAVID.p.value, y = InnateDB.p.value))+
  geom_point(aes(),alpha=0.4, size=4)

#zoom in

ggplot(ID50UP, aes(x = DAVID.p.value, y = InnateDB.p.value))+
  geom_point(aes(),alpha=0.4, size=4)+
  scale_x_continuous(lim=c(0,0.0001))+
  scale_y_continuous(lim=c(0,0.0001))



#How many terms in each day contain "mitotic"or "mitosis"?

dayList2<-c("4","7","14")

a<-function(Day){
  b<-ID50UP[ID50UP$Day==Day,]
sum(str_count(b$Term,"mitotic")|str_count(b$Term,"mitosis"))
}
c<-lapply(dayList2,FUN=a)

#there are 13 unique terms containing mitosis or mitotic

d<-unique(ID50UP$Term[str_detect(ID50UP$Term,
                              "mito(s|t)i(s|c)")==TRUE])


length(unique(ID50UP$Term))

################DAVID + Innate data (not1) for dose= 50 DOWN
overlapInnate50not1DOWN<-combinedInnate50 %>%
  filter(direction=="DOWN", Day!="1")%>%
  group_by(Pathway.Id)%>%
  filter(n()==3)%>%
  ungroup()

# how many go terms are common to all 3 days
overlapInnate50not1DOWN %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0)) %>%
   ungroup()  %>% 
   summarize(allSig = sum(allSig))
# says 7

#separating go terms out from the pre-existing david data
overlapDnot150DOWN <- subsetToOverlappingGoTerms(getAllDavid(50), "DOWN", 
   withDay1 = FALSE)
overlapDnot150DOWN<-getGO(overlapDnot150DOWN)

#merge the data sets
ID50DOWN<-merge(overlapInnate50not1DOWN,overlapDnot150DOWN,
              by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID50DOWN<-ID50DOWN%>%
  select(Day,Pathway.Id,Term,Pathway.p.value.corrected,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value.corrected<0.05)

names(ID50DOWN)[4:5]= c("InnateDB.p.value","DAVID.p.value")

ID50DOWN$Day<-factor(ID50DOWN$Day, levels = c("4","7","14"))

length(unique(ID50DOWN$Term))

#see DOI: 10.1111/hiv.12100 on gingival tissue tenofovir


length(unique(ID50DOWN$Term[str_detect(ID50DOWN$Term,
                              "mito(s|t)i(s|c)")==TRUE]))

################ Innate data for dose = 500
combinedInnate500 <- getAllInnateDB(500)

################DAVID + Innate data (not1) for dose= 500 UP
overlapInnate500not1UP<-combinedInnate500 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Pathway.Id)%>%
  filter(n()==3)%>%
  ungroup()

# how many go terms are common to all 3 days
overlapInnate500not1UP %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0)) %>%
   ungroup()  %>% 
   summarize(allSig = sum(allSig))
# says 47

overlapInnate500not1UP %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0),
      z = Pathway.Name[1]) %>%
   filter(allSig == 1)

#separating go terms out from the pre-existing david data
overlapDnot1500up <- subsetToOverlappingGoTerms(getAllDavid(500), "UP", 
   withDay1 = FALSE)
overlapDnot1500up<-getGO(overlapDnot1500up)

#merge the data sets
ID500UP<-merge(overlapInnate500not1UP,overlapDnot1500up,
                by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID500UP<-ID500UP%>%
  select(Day,Pathway.Id,Term,Pathway.p.value.corrected,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value.corrected<0.05)

names(ID500UP)[4:5]= c("InnateDB.p.value","DAVID.p.value")

ID500UP$Day<-factor(ID500UP$Day, levels = c("4","7","14"))

length(unique(ID500UP$Term))

length(unique(ID500UP$Term[str_detect(ID500UP$Term,
                                       "mito(s|t)i(s|c)")==TRUE]))



################DAVID + Innate data (not1) for dose= 500 down
overlapInnate500not1DOWN<-combinedInnate500 %>%
   filter(direction=="DOWN", Day!="1")%>%
   group_by(Pathway.Id)%>%
   filter(n()==3)%>%
   ungroup()

# how many go terms are common to all 3 days
overlapInnate500not1DOWN %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0)) %>%
   ungroup()  %>% 
   summarize(allSig = sum(allSig))
# says 10

overlapInnate500not1DOWN %>%
   group_by(Pathway.Id) %>%
   summarize(allSig = ifelse(all(Pathway.p.value.corrected < 0.05), 1, 0),
      z = Pathway.Name[1]) %>%
   filter(allSig == 1)
   
   


############ biomaRt#############################
######## Finding # of genes assoc with GO ids
require(illuminaHumanv4.db)
require(biomaRt)
require(GO.db)
library(topGO)
library(GOstats)
library(annotate)
library(genefilter)

# I want to make a column in a df of overlap data with the
#number of genes in each go id.

#First make a column for the go terms using getGO
#NOTE: if there is no GO term, there will be an NA there
allD50<-getGO(allD50)

#make a list of unique go Ids from allD50
allD50GOlist<-as.list(unique(allD50$Pathway.Id))

#mapping of Go ids to entrez ids for my universe
#Using topGO functions and methods
GOID2Gene<-annFUN.org(c("CC","BP","MF"),
                      feasibleGenes = NULL,
                                  mapping="org.Hs.eg.db",
                      ID = "entrez")

#geneNames are the GO ids associated with entrez Ids
#they are the names of the elemnts in the GOID2Gene list
geneNames<-names(GOID2Gene)

#filter the elements in the universe that overlap
#with list of GO ids from allD50.
geneList<-GOID2Gene[geneNames %in% allD50GOlist]

#here is a named list where the names are GO ids from allD50
#and the element is the number of associated genes
GOidLength<-lapply(geneList,length)

#convert the GOidLength list into a dataframe
#used dplyr's as_data_frame so the GO ids aren't turned into
#rownames
GOidLength<-as_data_frame(GOidLength)

#now it is wide and I want it to be long so I will melt
GOidLength<-melt(GOidLength)

#cleanup so it will merge well
colnames(GOidLength)<-c("Pathway.Id","GenesInGOid")
GOidLength$Pathway.Id<-as.character(GOidLength$Pathway.Id)
GOidLength$GenesInGOid<-as.numeric(GOidLength$GenesInGOid)

#merge it with the allD50 dataframe
allD50<-merge(GOidLength,allD50, by = "Pathway.Id")

### NOTE: This no longer contains data from non-GO databases
### like SPIR and etc.










