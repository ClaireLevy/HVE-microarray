
require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
require(pander)
source("getAllDavid.R")

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

overlapD50up<-allD50%>%
  filter(direction=="UP")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


#if you do summarize(n())where filter is  it shows you the term and 
#how many occurances there were in overlapTerms.50.UP. We only want the terms where
#there there 4 occurances (1 per day we looked at)


overlapD50upshort<-overlapD50up%>%
  select(Day, Category, Term, Benjamini)%>%
  filter(Benjamini<0.05)


#a function to cast the data into an easy to read format
castFunction<-function(df){
  dcast(df, Category + Term~Day,
        value.var = "Benjamini")
}
overlapD50upshort<-castFunction(overlapD50upshort)
#get rid of NAs
overlapD50upshort<-na.omit(overlapD50upshort)

#note that all the day 1 values have dropped out becasue
#they were all >0.05 but the terms still reflect overlaps 
#in all 4 days.

write.csv(overlapD50upshort,
          "overlap.DAVID.50.UP.csv", row.names=FALSE)


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

overlapDnot150up<-allD50 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDnot150upshort<-overlapDnot150up %>%
         select(Day,Category,Term,Benjamini)

overlapDnot150upshort<-castFunction(overlapDnot150upshort)

# function that takes a data frame with the following columns, in order, 
# "Category", "Term", "4", "7", "14" 
# and returns the rows where p-values are < 0.05 on all three days
allSig <- function(data) {
   # check data for validity
   if (any(colnames(data) != c("Category", "Term", "4", "7", "14"))) {
      cols <- paste(colnames(data), collapse = ", ")
      stop(paste("allSig expects the column names of data to be: ", 
         "'Category, Term, 4, 7, 14', but the column names of data are '",
         cols, "' so this might be a mistake.", sep = ""))
   }
      
   # rename columns because dplyr gets mad at columns named numbers
   colnames(data) <- c("Category", "Term", "d4", "d7", "d14")
   
   # remove rows where one or more p-value is > 0.05
   data %>%
      filter(d4 < 0.05, d7 < 0.05, d14 < 0.05)   
}

allSigNot150Up <- allSig(overlapDnot150upshort)

write.csv(allSigNot150Up,
          file = "overlap.DAVID.not1.50.UP.csv", row.names = FALSE)


###############DAVID data for dose = 50 DOWN ############
# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

overlapD50DOWN<-allD50 %>%
  filter(direction=="DOWN")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


ggplot(data=overlapD50DOWN, aes())+
  geom_point(aes(x = Day , y = Benjamini),
             position=position_jitter(w=0.15),size=3)+
  geom_hline(y=0.05)+
  ggtitle("P-values for overlapping GO terms \n\ in DOWN-regulated genes,dose=50")



overlapD50DOWNshort<-overlapD50DOWN%>%
  select(Day, Category, Term, Benjamini)%>%
  filter(Benjamini<0.05)
  
  
#note that all the day 1 values have dropped out becasue
#they were all >0.05 but the terms still reflect overlaps 
#in all 4 days.

overlapD50DOWNshort<-na.omit(castFunction(overlapD50DOWNshort))


write.csv(overlapD50DOWNshort,
          file = "overlap.DAVID.50.DOWN.csv", row.names=FALSE)



##############DAVID data for dose = 50 DOWN not incl day 1##############
# read in the files from DAVID for concentration == 50
allD50 <- getAllDavid(50)

overlapDnot150DOWN<-allD50 %>%
  filter(direction=="DOWN", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDnot150DOWNshort<-overlapDnot150DOWN %>%
  select(Day,Category,Term,Benjamini)

overlapDnot150DOWNshort<-castFunction(overlapDnot150DOWNshort)

allSigNot150Down <- allSig(overlapDnot150DOWNshort)


write.csv(allSigNot150Down,
          file = "overlap.D.not1.50.DOWN.csv", row.names = FALSE)


############D data for dose = 500 UP #############
# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)

overlapD500up<-allD500%>%
  filter(direction=="UP")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


overlapD500upshort<-overlapD500up%>%
  select(Day, Category, Term, Benjamini)%>%
  filter(Benjamini<0.05)


#use castFunction
overlapD500upshort<-castFunction(overlapD500upshort)

overlapD500upshort<-na.omit(overlapD500upshort)

write.csv(overlapD500upshort,
          "overlap.DAVID.500.UP.csv",row.names=FALSE)


###############DAVID data for dose = 500 UP excluding day 1##############
# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)

overlapDnot1500up<-allD500 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDnot1500upshort<-overlapDnot1500up %>%
  select(Day,Category,Term,Benjamini)

overlapDnot1500upshort<-castFunction(overlapDnot1500upshort)

allSigNot1500Up <- allSig(overlapDnot1500upshort)

write.csv(allSigNot1500Up,
          file = "overlap.D.not1.500.UP.csv", row.names = FALSE)


###############DAVID data for dose = 500 DOWN ###########
##################DAVID DATA NEEDS TO BE CORRECTED##########

# read in the files from DAVID for concentration == 500
allD500 <- getAllDavid(500)


#####excluding day 1##############



########################Innate DB data#######################

#########Innate data dose = 50###############################

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50")

innateFiles50<-list.files("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50",
                          pattern="^Innate")

innateData50<-lapply(innateFiles50,read.csv, sep="\t",quote="")

names(innateData50) <- stringr::str_replace(innateFiles50, pattern = ".csv", replacement = "")


innateData50<-Map(cbind,innateData50,Day=dayList,direction=directionList)

combinedInnate50<-rbind_all(innateData50)



##function to separate out the GO ids in the DAVID data
##so they can be compared with innate data
getGO<-function(df){
  mutate(df,
         Pathway.Id= ifelse(str_detect(df$Term,"GO:")==TRUE,
                            substr(df$Term,1,10),NA))
}

#apply the getGO function over them so they all have GO ids
Ddata50<-lapply(Ddata50,getGO)

combinedD50<-rbind_all(Ddata50)

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

#separating go terms out from the pre-existing david data
overlapDnot150up<-getGO(overlapDnot150up)

#merge the data sets
ID50UP<-merge(overlapInnate50not1UP,overlapDnot150up,
         by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID50UP<-ID50UP%>%
  select(Day,Pathway.Id,Term,Pathway.p.value..corrected.,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value..corrected.<0.05)

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

#separating go terms out from the pre-existing david data
overlapDnot150DOWN<-getGO(overlapDnot150DOWN)

#merge the data sets
ID50DOWN<-merge(overlapInnate50not1DOWN,overlapDnot150DOWN,
              by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID50DOWN<-ID50DOWN%>%
  select(Day,Pathway.Id,Term,Pathway.p.value..corrected.,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value..corrected.<0.05)

names(ID50DOWN)[4:5]= c("InnateDB.p.value","DAVID.p.value")

ID50DOWN$Day<-factor(ID50DOWN$Day, levels = c("4","7","14"))

length(unique(ID50DOWN$Term))

#see DOI: 10.1111/hiv.12100 on gingival tissue tenofovir


length(unique(ID50DOWN$Term[str_detect(ID50DOWN$Term,
                              "mito(s|t)i(s|c)")==TRUE]))

################ Innate data for dose = 500

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 500")

innateFiles500<-list.files("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 500",
                          pattern="^Innate")

innateData500<-lapply(innateFiles500,read.csv, sep="\t",quote="")

names(innateData500) <- stringr::str_replace(innateFiles500, pattern = ".csv", replacement = "")


innateData500<-Map(cbind,innateData500,Day=dayList,direction=directionList)

combinedInnate500<-rbind_all(innateData500)

################DAVID + Innate data (not1) for dose= 500 UP
overlapInnate500not1UP<-combinedInnate500 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Pathway.Id)%>%
  filter(n()==3)%>%
  ungroup()

#separating go terms out from the pre-existing david data
overlapDnot1500up<-getGO(overlapDnot1500up)

#merge the data sets
ID500UP<-merge(overlapInnate500not1UP,overlapDnot1500up,
                by=c("Day","Pathway.Id"))



#select only <0.05 p vals
ID500UP<-ID500UP%>%
  select(Day,Pathway.Id,Term,Pathway.p.value..corrected.,
         Benjamini)%>%
  filter(Benjamini<0.05 & Pathway.p.value..corrected.<0.05)

names(ID500UP)[4:5]= c("InnateDB.p.value","DAVID.p.value")

ID500UP$Day<-factor(ID500UP$Day, levels = c("4","7","14"))

length(unique(ID500UP$Term))

length(unique(ID500UP$Term[str_detect(ID500UP$Term,
                                       "mito(s|t)i(s|c)")==TRUE]))
############ biomaRt#############################
#Goal: get GO ids and terms from bioMaRt from list of entrez ids

require(illuminaHumanv4.db)
require(biomaRt)

UP<-longForm%>%
  filter(Direction == "UP", ENTREZ_GENE_ID)
#make a vector of the ids to annotate
UPentrez<-UP$ENTREZ_GENE_ID

#set the mart you want to use (see listMarts())
#and choose the dataset to use see listDatasets()

ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#attributes: the values you want
#filters: the data biomaRt should look at and retrieve attributes for   

#the attribute for go term name is "name_1066"


#This is a list of the entrez ids and the associated go terms and IDs
goids<-getBM(attributes=c("entrezgene","go_id","name_1006"),
            filters = "entrezgene",values=UPentrez,
            mart = ensembl)


