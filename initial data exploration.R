
require(dplyr)
require(reshape2)
require(stringr)
dataSummary<-read.csv("ex_vivo_HVE_Tenofovir_Summary.csv")
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


##This function doesn't work but I want it to :(
#look more at ?interp
isolateDEG<-function("colname","direction"){
summaryOnly%>%
  select("Probe.ID","TargetID","ENTREZ_GENE_ID","colname")%>%
  filter(colname=="direction")
}



  

###all the DOWNs

d14.50.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d14.50)%>%
  filter(d14.50=="DOWN")

d7.50.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d7.50)%>%
  filter(d7.50=="DOWN")

d1.50.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d1.50)%>%
  filter(d1.50=="DOWN")

d14.500.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d14.500)%>%
  filter(d14.500=="DOWN")

d7.500.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d7.500)%>%
  filter(d7.500=="DOWN")

d1.500.down<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d1.500)%>%
  filter(d1.500=="DOWN")



#all the UPs

d14.50.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d14.50)%>%
  filter(d14.50=="UP")

d7.50.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d7.50)%>%
  filter(d7.50=="UP")

d1.50.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d1.50)%>%
  filter(d1.50=="UP")

d14.500.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d14.500)%>%
  filter(d14.500=="UP")

d7.500.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d7.500)%>%
  filter(d7.500=="UP")

d1.500.up<-summaryOnly%>%
  select(Probe.ID,TargetID,ENTREZ_GENE_ID,d1.500)%>%
  filter(d1.500=="UP")




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



ggplot(meltedupDownCount,aes(x=day.dose,y=Count))+
  geom_bar(aes(fill=Direction), position="dodge",stat="identity")+ylab("Number of genes")+
  geom_hline(y=7832)+
  ggtitle("Number of upregulated, downregulated and non-significant genes\n\ per time point and dose\n\ 
          (line represents total genes analyzed for diff exp)")