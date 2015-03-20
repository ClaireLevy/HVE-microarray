
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

###############DAVID data for dose = 50 UP ######################################
#NOTE: on DAVID site I used the illumina HT 12 v3 background 
# and  in the functional annotation chart,under options
#checked "Fisher exact", then reran with options and saved.

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray")

#extract the david data from the folder where I saved it

#make a list of the files I want
DAVIDlist50<-list.files("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50",
                        pattern = "^DAVID")

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50")

#apply read.csv over the list
DAVIDdata50<-lapply(DAVIDlist50,read.csv)

#give the list elements names

names(DAVIDdata50)<- str_replace(DAVIDlist50, pattern = ".csv", replacement = "")

#check the order by looking at names(DAVIDdata)
#make a list for day and direction 
dayList<-list(1,1,14,14,4,4,7,7)
directionList<-list("DOWN","UP")

# add columns to the dfs for day and direction
DAVIDdata50<-Map(cbind,DAVIDdata50,Day=dayList)
DAVIDdata50<-Map(cbind,DAVIDdata50,direction=directionList)

#make the list of dfs into one df

allDAVID50<-dplyr::rbind_all(DAVIDdata50)

overlapDAVID50up<-allDAVID50%>%
  filter(direction=="UP")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


#if you do summarize(n())where filter is  it shows you the term and 
#how many occurances there were in overlapTerms.50.UP. We only want the terms where
#there there 4 occurances (1 per day we looked at)


overlapDAVID50upshort<-overlapDAVID50up%>%
  select(Day, Category, Term, Benjamini)


#a function to cast the data into an easy to read format
castFunction<-function(df){
  dcast(df, Category + Term~Day,
        value.var = "Benjamini")
}
overlapDAVID50upshort<-castFunction(overlapDAVID50upshort)


write.csv(overlapDAVID50upshort, "overlap.DAVID.50.UP.csv")


###plot plot plot

#Pvalues plot
plot150UP<-ggplot(data=overlapDAVID50up, aes())+
  geom_point(aes(x = Day , y = Benjamini),
             position=position_jitter(w=0.15),
             size=3, 
             alpha = 0.7)+
  geom_hline(y=0.05)+
  ggtitle("P-values for overlapping GO terms \n\ in Up-regulated genes,dose=50")
  

# This will save a 400x400 file at 100 ppi
ggsave("plot150UP.png", width=4, height=4, dpi=100)


#########DAVID data for dose = 50 UP not incl day 1###########\
overlapDAVIDnot150up<-allDAVID50 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDAVIDnot150upshort<-overlapDAVIDnot150up %>%
         select(Day,Category,Term,Benjamini)

overlapDAVIDnot150upshort<-castFunction(overlapDAVIDnot150upshort)

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

allSigNot150Up <- allSig(overlapDAVIDnot150upshort)

write.csv(allSigNot150Up,
          file = "overlap.DAVID.not1.50.UP.csv", row.names = FALSE)

#pvalues for those overlapping terms, more low ones for d7
ggplot(data=overlap.DAVID.not1.50.UP, aes())+
  geom_point(aes(x = Day , y = Benjamini),
             position=position_jitter(w=0.15),size=3,
             alpha=0.2)+
  geom_hline(y=0.05)+
  ggtitle("Uncorrected Pvalues of overlapping upreg terms in day 4,7,& 14
          \n\ dose = 50")

###############DAVID data for dose = 50 DOWN ############
#extracting data from the long form summary
extract("1","50","DOWN")
extract("4","50","DOWN")
extract("7","50","DOWN")
extract("14","50","DOWN")


overlapDAVID50DOWN<-allDAVID50 %>%
  filter(direction=="DOWN")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()


ggplot(data=overlapDAVID50DOWN, aes())+
  geom_point(aes(x = Day , y = Benjamini),
             position=position_jitter(w=0.15),size=3)+
  geom_hline(y=0.05)+
  ggtitle("P-values for overlapping GO terms \n\ in DOWN-regulated genes,dose=50")



overlapDAVID50DOWN<-castFunction(select(overlapDAVID50DOWN,
                                  Day,Category,Term,Benjamini))



write.csv(overlapDAVID50DOWN,
          file = "overlap.DAVID.50.DOWN.csv")
##############DAVID data for dose = 50 DOWN not incl day 1##############
overlapDAVIDnot150DOWN<-allDAVID50 %>%
  filter(direction=="DOWN", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDAVIDnot150DOWNshort<-overlapDAVIDnot150DOWN %>%
  select(Day,Category,Term,Benjamini)

overlapDAVIDnot150DOWNshort<-castFunction(overlapDAVIDnot150DOWNshort)

allSigNot150Down <- allSig(overlapDAVIDnot150DOWNshort)


write.csv(allSigNot150Down,
          file = "overlap.DAVID.not1.50.DOWN.csv", row.names = FALSE)


############DAVID data for dose = 500 UP #############
DAVIDlist500<-list.files("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 500",
                        pattern = "^DAVID")

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 500")

#apply read.csv over the list
DAVIDdata500<-lapply(DAVIDlist500,read.csv)

#give the list elements names

names(DAVIDdata500)<- str_replace(DAVIDlist500, pattern = ".csv", replacement = "")


# add columns to the dfs for day and direction
#using dayList and directionList from above
DAVIDdata500<-Map(cbind,DAVIDdata500,Day=dayList)
DAVIDdata500<-Map(cbind,DAVIDdata500,direction=directionList)

#make the list of dfs into one df

allDAVID500<-dplyr::rbind_all(DAVIDdata500)

overlapDAVID500up<-allDAVID500%>%
  filter(direction=="UP")%>%
  group_by(Term)%>%
  filter(n()==4)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDAVID500upshort<-overlapDAVID500up%>%
  select(Day, Category, Term, Benjamini)


#use castFunction
overlapDAVID500upshort<-castFunction(overlapDAVID500upshort)


write.csv(overlapDAVID500upshort, "overlap.DAVID.500.UP.csv")


###############DAVID data for dose = 500 UP excluding day 1##############
overlapDAVIDnot1500up<-allDAVID500 %>%
  filter(direction=="UP", Day!="1")%>%
  group_by(Term)%>%
  filter(n()==3)%>% # because there are 4 days
  arrange(Day,Term,PValue)%>%
  ungroup()



overlapDAVIDnot1500upshort<-overlapDAVIDnot1500up %>%
  select(Day,Category,Term,Benjamini)

overlapDAVIDnot1500upshort<-castFunction(overlapDAVIDnot1500upshort)

allSigNot1500Up <- allSig(overlapDAVIDnot1500upshort)

write.csv(allSigNot1500Up,
          file = "overlap.DAVID.not1.500.UP.csv", row.names = FALSE)


###############DAVID data for dose = 500 DOWN ###########
##################DAVID DATA NEEDS TO BE CORRECTED##########
extract("1","500","DOWN")
extract("4","500","DOWN")
extract("7","500","DOWN")
extract("14","500","DOWN")


#####excluding day 1##############



#############################################################
########################Innate DB data#######################

#################### I tried this but it didn't work :( 

#Innate.50.UPfiles<-file.path("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50",
                             c("InnateDB-1.50.up.txt","InnateDB-4.50.up.txt",
                "InnateDB-7.50.up.txt","InnateDB-14.50.up.txt"))


#Innate.50.UP<-lapply(Innate.50.UPfiles,read.table, sep="\t")
########################################################

getInnate<-function(day,concentration,direction){
  read.table(file = paste("InnateDB",day,concentration,direction,"txt", sep="."),
           sep = "\t",header=TRUE, quote="")
}
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes/dose = 50")

Innate.1.50.UP<-getInnate("1","50","UP")

Innate.4.50.UP<-getInnate("4","50","UP")

Innate.7.50.UP<-getInnate("7","50","UP")
 
Innate.14.50.UP<-getInnate("14","50","UP")
                      

##function to separate out the GO ids in the DAVID data
##so they can be compared with innate data
getGO<-function(df){
  mutate(df,
         Pathway.Id= ifelse(str_detect(df$Term,"GO:")==TRUE,
                            substr(df$Term,1,10),NA))
}

#make a list of the david data frames
DAVID.50.UPdfs<-list(DAVID.1.50.UP,DAVID.4.50.UP,DAVID.7.50.UP,
              DAVID.14.50.UP)

#apply the getGO function over them so they all have GO ids
DAVID.50.UPdfs<-lapply(DAVID.50.UPdfs,getGO)



#which GO ids are in both innate and david data for day 1 dose = 50?

#USe the InnateDAVID function to merge

setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/differentiallyExpressedGenes")
source("InnateDAVIDoverlap.R")


#(Subset from the DAVID dfs list to refer to the one
#you want)

ID.1.50.UP<-InnateDAVID(DAVID.50.UPdfs[1],Innate.1.50.UP,50)

ID.4.50.UP<-InnateDAVID(DAVID.50.UPdfs[2],Innate.4.50.UP,50)

ID.7.50.UP<-InnateDAVID(DAVID.50.UPdfs[3],Innate.7.50.UP,50)

ID.14.50.UP<-InnateDAVID(DAVID.50.UPdfs[4],Innate.14.50.UP,50)         
         
#combine all of the merged (overlapping GO term) data          
allID.50.UP<-rbind(ID.1.50.UP,ID.4.50.UP,
                   ID.7.50.UP,ID.14.50.UP)

names(allID.50.UP)[3:4]<-c("InnateDB.p.value","DAVID.p.value")


#calculate the abs value of the difference between p values
#and arrange by day and that difference
#Terms are the top of the list for each day have 
#most similar p values. Both are Benj corrected
#I THINK DAVID data is a correction of fisher's exact
#Innate is a correction of Hypergeometric
  
allID.50.UP<-allID.50.UP %>%
  mutate(abs.difference = abs(InnateDB.p.value-DAVID.p.value))%>%
  arrange(Day,abs.difference)

allID.50.UP$Day<-factor(allID.50.UP$Day, levels=c("1","4","7","14"))

#InnateDB and DAVID p values
ggplot(allID.50.UP, aes(x = Day, y = value, color = variable))+
  geom_point(aes(y=InnateDB.p.value,
                 col="InnateDB.p.value"),
             position=position_jitter(w=0.15),alpha=0.5)+
  geom_point(aes(y=DAVID.p.value, col="DAVID.p.value"),
             position=position_jitter(w=0.15),alpha=0.5)+
  labs(y="Benjamini")+
  theme(legend.title=element_blank())+
  ggtitle("Overlapping GO terms from InnateDB and DAVID ORA \n\
up-regulated genes, dose = 50")

#InnateDB VS DAVID pvalues, limited 0-0.05
ggplot(allID.50.UP, aes(x = DAVID.p.value, y = InnateDB.p.value))+
  geom_point(aes())+
  scale_x_continuous(limits=c(0,0.05))+
  scale_y_continuous(limits= c(0,0.05))


write.csv(allID.50.UP,"allID.50.UP.csv")




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

# Remove genes with no Entrez

entrezIds<-mget()

#This has the potential to minimize all the txt and csv files..
