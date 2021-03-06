setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray")
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


# old plotting code for ORA
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

###############################################################################
# start of DAVID / InnateDB ORA 
###############################################################################
source("getAllDavid.R")
source("getAllInnateDB.R")
source('subsetToOverlappingGoTerms.R')
#source("write-files-for-DAVID-and-InnateDB.R") # rewrites files, don't need each time
source("analyzeAndWriteGO.R")
#################DAVID/INNATE ORA OVERALL #####################
# The overall overrepresentation analysis involves:
# 1. Read in LMF's DEG spreadsheet
# 2. Convert it to .txt and .csv files with the DEGs for each
#     day, concentration of tenofovir, and direction of change
# 3. Enter those files into InnateDB and DAVID
# 4. Read in the files from InnateDB and DAVID containing the
#     overrepresented GO and other terms
# 5. Find terms that are enriched on multiple days
# 6. Write out summary files with the results.

#################DAVID/INNATE ORA INPUT #######################

# Steps 1-2
# To write .txt and .csv files for input into DAVID/InnateDB,
# use write-files-for-DAVID-and-InnateDB.R

# This code converts LMF's DEG spreadsheet into longform df
# and converts that to .txt and .csv files to enter into
# InnateDB and DAVID. 
#WHAT TO PUT IN: nothing
#WHAT YOU GET: longform df (includes DEG's only)
#.txt and .csv files in dose = 500 or 50 folders


#################DAVID/INNATE WEBSITE #######################

# Step 3
# See getAllDavid for details of website use

###############GET DAVID RESULTS AS DF #####################

# Step 4
#WHAT TO PUT IN: desired concentration (50 or 500)
#WHAT YOU GET: df with all DAVID data for that concentration


allD50 <- getAllDavid(50)
allD500 <- getAllDavid(500)

###############GET INNATEDB RESULTS AS DF #####################

# Step 4
#WHAT TO PUT IN: desired concentration (50 or 500)
#WHAT YOU GET: df with all InnateDB data for that concentration, with columns
#              added to match DAVID template


allInnate50<-getAllInnateDB(50)
allInnate500<-getAllInnateDB(500)

########## FIND GO TERMS THAT OVERLAP ON MULTIPLE DAYS######

# Steps 5-6
#FUNCTION: analyzeAndWriteGO
#WHAT TO PUT IN: 
#        df:            allD50(0) or allInnate50(0)
#        concentration: 50 or 500
#        filename:      e.g. "overlap.DAVID.not1.50.DOWN.csv"
#        direction:     "UP" or "DOWN"
#WHAT YOU GET: a CSV file containing the terms that are enriched on days 4, 7, 
#              and 14 and p < 0.05. The file has the given name and is in the 
#              appropriate folder 

analyzeAndWriteGO(allD50, 50, "overlap.DAVID.not1.50.DOWN.csv", "DOWN")
analyzeAndWriteGO(allD50, 50, "overlap.DAVID.not1.50.UP.csv", "UP")
analyzeAndWriteGO(allD500, 500, "overlap.DAVID.not1.500.DOWN.csv", "DOWN")
analyzeAndWriteGO(allD500, 500, "overlap.DAVID.not1.500.UP.csv", "UP")

analyzeAndWriteGO(allInnate50, 50, "overlap.Innate.not1.50.DOWN.csv", "DOWN")
analyzeAndWriteGO(allInnate50, 50, "overlap.Innate.not1.50.UP.csv", "UP")
analyzeAndWriteGO(allInnate500, 500, "overlap.Innate.not1.500.DOWN.csv", "DOWN")
analyzeAndWriteGO(allInnate500, 500, "overlap.Innate.not1.500.UP.csv", "UP")

# Step 5-6
#FUNCTION: analyzeAndWriteGOWithDay1
#WHAT TO PUT IN: 
#        df:            allD50(0) or allInnate50(0) (from getAllX())
#        concentration: 50 or 500
#        filename:      e.g. "overlap.DAVID.not1.50.DOWN.csv"
#        direction:     "UP" or "DOWN"
#        day1Sig:       whether day 1 p-value must be significant (TRUE/FALSE)
#WHAT YOU GET: a CSV file containing the terms that are enriched on days 1, 4, 7, 
#              and 14 and p < 0.05 for 4, 7, 14, and optionally also for d1.
#              The file has the given name and is in the appropriate folder 


analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.DOWN.csv", "DOWN")
analyzeAndWriteGOWithDay1(allD50, 50, "overlap.DAVID.50.UP.csv", "UP")
analyzeAndWriteGOWithDay1(allD500, 500, "overlap.DAVID.500.DOWN.csv", "DOWN", TRUE)
analyzeAndWriteGOWithDay1(allD500, 500, "overlap.DAVID.500.UP.csv", "UP", TRUE)

analyzeAndWriteGOWithDay1(allInnate50, 50, "overlap.Innate.50.DOWN.csv", "DOWN")
analyzeAndWriteGOWithDay1(allInnate50, 50, "overlap.Innate.50.UP.csv", "UP")
analyzeAndWriteGOWithDay1(allInnate500, 500, "overlap.Innate.500.DOWN.csv", "DOWN", TRUE)
analyzeAndWriteGOWithDay1(allInnate500, 500, "overlap.Innate.500.UP.csv", "UP", TRUE)

# Step 5 (optional)
#FUNCTION: subsetToOverlappingGoTerms
#WHY? This function is automatically used by the analyzeAndWriteGO functions, so
#     it's not necessary to use it manually to generate the summary files. 
#     However, for other purposes, it can be convenient to use to get the terms.
#WHAT TO PUT IN: 
#           df:         A df of innate or david data (use getAllX)
#           direction:  "UP" or "DOWN"
#           withDay1:   designation of T or F for incl day1
#WHAT YOU GET: A df containing the GO (and for DAVID, other) terms that 
#           overlapped between days 4, 7, and 14, and optionally day 1. It 
#           includes all the terms that overlapped, not just those that had 
#           significant p-values.  

source("subsetToOverlappingGoTerms.R")

###############################################################################
# end of DAVID / InnateDB ORA 
###############################################################################


####DAVID + Innate data (not1) for dose = 50 UP
#DAVID overlaps
overlapDnot150up <- subsetToOverlappingGoTerms(getAllDavid(50),
                                               "UP", 
                                               withDay1 = FALSE)
#Innate overlaps

overlapInnatenot150up <- subsetToOverlappingGoTerms(getAllInnateDB(50),
                                               "UP", 
                                               withDay1 = FALSE)

















