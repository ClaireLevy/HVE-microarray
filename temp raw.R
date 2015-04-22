require(dplyr)
library(lumi)
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData")
RAW<-'FinalReport_032615_exvivo_TNF_HVE.txt'
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

RAWlumi.T = lumiT(RAW.lumi)
RAWlumi.N = lumiN(RAWlumi.T, method = "rsn")
RAWlumi.N.Q = lumiQ(RAWlumi.N, detectionTh = 0.05)

###################Filter undetected/low SD genes ##############
library(limma)
require(genefilter)
dataMatrix<-exprs(RAWlumi.N.Q)

#remove unexpressed and and un-annotated probes
presentCount<-detectionCall(RAWlumi.N.Q,Th=0.05, type = "probe")
selDataMatrix<-dataMatrix[presentCount>0,]
nrow(dataMatrix)
#> [1] 47323
nrow(selDataMatrix)
#> [1] 32921

#remove probes with low row SD
selDataMatrixSds<-rowSds(selDataMatrix)
selDataMatrixFiltered<-selDataMatrix[selDataMatrixSds>=0.2,]
nrow(selDataMatrixFiltered)
#> [1] 10421

hist(dataMatrix,breaks=50)
hist(selDataMatrix, breaks=50)
hist(selDataMatrixFiltered, breaks=50)

################ FIT MODEL ################################
dose500day14Data<-selDataMatrixFiltered[,c("HVE_A4", "HVE_A8","HVE_A12",
                                   "HVE_C4", "HVE_C8","HVE_C12")]

design<-model.matrix(~0+
      factor(c(0, 0, 0, 500, 500, 500), levels = c("0", "500")))
colnames(design)<-c("Control","Treatment")
design
#>   Control Treatment
#> 1       1         0
#> 2       1         0
#> 3       1         0
#> 4       0         1
#> 5       0         1
#> 6       0         1

#fit the model
fit<-lmFit(dose500day14Data,design)

# fit contrast of interest
cont.matrix<-makeContrasts(TreatmentvsControl=Treatment-Control,
                           levels=design)
cont.matrix
#> Levels      TreatmentvsControl
#>   Control                   -1
#>   Treatment                  1

fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
topTable(fit2,adjust="BH")
#> adjusted P value < 0.05 for 0 probes

#non-contrasts matrix?
fit<-eBayes(fit)
topTable(fit,adjust="BH")
#> adjusted P value < 0.05 for ALL probes