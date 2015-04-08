#Working with raw illumina data
#started 26March15


require(dplyr)
library(lumi)
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData")
#LMF sent FinalReport_exvivoTNF_HVE_032615.txt which is raw
#data with control PROBE profile data rather than control
#GENE profile data, which we tried before and couldn't get lumi
#to read.

#she also sent her code for reading in and normalising (below)

### LMF:load data into memory 
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
#why use lumiHumanAll.db instead of illuminaHumanv4.db?

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

###############################################################################
# explore sample relations
###############################################################################

# colors 
red <- "#E41A1C"
blue <- "#377EB8"
green <- "#4DAF4A"
purple <- "#984EA3"

# get key to experiment IDs
key <- read.csv("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData/SampleKey.csv",
   header = T)
key$MicroarrayId <- paste("HVE_", key$MicroarrayId, sep = "")

# sort the "key" data frame according to the order that the samples are in in 
# the lumi file for convenience when providing the color arguments
sortKey <- function() {
   key[match(sampleNames(RAWlumi.N.Q), key$MicroarrayId), ]
}

# convenience plot function
p <- function(title) {
   lumi::plotSampleRelation(RAWlumi.N.Q, method="mds", color = sortKey()$Color, 
      main = title)
}

# convenience plot function (limma)
p2 <- function(title) {
   limma::plotMDS(RAWlumi.N.Q, pch = 16, col = sortKey()$Color, 
      main = title, top = 3000)
   # limma bases the distance on "top" number of genes
   # I chose 3000 because it's roughly how many were differentially expressed
   
   # for some reason it's the mirror image of the lumi plots
}

# by TFV concentration
key$Color <- ifelse(key$Concentration == 0, red, 
   ifelse(key$Concentration == 50, green, blue))
p("Colors are messed up damn you lumi Blue = 500 uM tenofovir, green = 50, red = 0")
p2("Blue = 500 uM tenofovir, green = 50, red = 0")

# by cell line
key$Color <- ifelse(key$CellLine == "HVE1", red, 
   ifelse(key$CellLine == "HVE2", green, blue))
p("Blue = HVE3, green = HVE2, red = HVE1")
p2("Blue = HVE3, green = HVE2, red = HVE1")

# by day
key$Color <- ifelse(key$Day == 1, red, 
   ifelse(key$Day == 4, green, 
      ifelse(key$Day == 7, blue, purple)))
p("Red = day 1, green = 4, blue = 7, purple = 14")
p2("Red = day 1, green = 4, blue = 7, purple = 14")

###################Continuing with lumi script, identify DEGs ##############
#Should I do the "inverse VST transform to the raw scale"??
#It seems like it is more for data with low expression values...

library(limma)
dataMatrix<-exprs(RAWlumi.N.Q)
head(dataMatrix)

#remove unexpressed and and un-annotated genes
#detectionCall "estimates the percentage of expressed genes
#of each sample" but I don't know what that means.
#I think this shows which probes have detection pvalues <0.01

presentCount<-detectionCall(RAWlumi.N.Q,Th=0.05, type = "probe")

#then this takes only those that have >0 occurances of pval<0.05
selDataMatrix<-dataMatrix[presentCount>0,]

#select columns with day 14 data according to exp design
#excel spreadsheet
dose500day14Data<-selDataMatrix[,c("HVE_A4", "HVE_A8","HVE_A12",
                                   "HVE_C4", "HVE_C8","HVE_C12")]

#making a design matrix based on limma vignette pg 41
#using the sample key for our experiment to get sample annotations

SampleKey<-read.csv("SampleKey.csv")

dose500Design<-SampleKey%>%
  dplyr::filter(Day==14, Concentration!=50)%>%
  dplyr::select(CellLine, Concentration)

Group<-factor(dose500Design$Concentration, levels=c("0","500"))
design<-model.matrix(~0+Group)#0 means no intercept I guess
colnames(design)<-c("Control","Treatment")


#fit the model
fit<-lmFit(dose500day14Data,design)

#Given a linear model fit to microarray data, compute estimated coefficients
#and standard errors for a given set of contrasts.
cont.matrix<-makeContrasts(TreatmentvsControl=Treatment-Control,
                           levels=design)


#Given a microarray linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical
#Bayes moderation of the standard errors towards a common value.

#fit model to contrasts matrix
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)

topTable(fit2,adjust="BH")

#non-contrasts matrix?
fit<-eBayes(fit)

topTable(fit,adjust="BH")


#get significant geen list with FDR adj p values<0.01
#lumi seems to skip the "fit2" step (no contrasts?)
require(lumiHumanAll.db)
require(annotate)
probeList<-rownames(selDataMatrix)
geneSymbol<-getSYMBOL(probeList,"illuminaHumanv4.db")
#lumi writes this whole thing as a function so there is sapply
#before lookUp
geneName<-lookUp(probeList,"illuminaHumanv4.db","GENENAME")

#make a df of probes, gene symbols and genenames
fit2$genes<-data.frame(ID=probeList, geneSymbol=geneSymbol,
                       geneName=geneName)

#print the top 10 genes
topTable(fit2,adjust="fdr",
               number=10)

p.adj<-p.adjust(fit$p.value[,2])

sigGene.adj<-probeList[p.adj<0.01]
