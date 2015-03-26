#Working with raw illumina data
#started 26March15



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

