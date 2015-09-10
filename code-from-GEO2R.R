library(lumi)
library(dplyr)
setwd("J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\microarrayData")

RAW<-'FinalReport_032615_exvivo_TNF_HVE.txt'
rawData <- lumiR(RAW,
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

transformed <- lumiT(rawData)
normalized <- lumiN(transformed, method = "rsn")
gset <- lumiQ(normalized, detectionTh = 0.05)

fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
# Sean: This code selects just the day 14, 500 uM samples. It depends on the
#        sample columns being in the order that they are on GEO/the order in
#        the file initially. The filtering code we wrote reorders them so this
#        code wouldn't work with the reordered data frame and would need to be
#        changed. 
sml <- c("X","X","X","X","X","G0","X","G0","X","G1","X","X","X","X","G0","X","X","X","X","X","X","X","X","X","G1","G1","X","X","X","X","X","X","X","X","X","X");

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# Sean: This code seems to be doing some sort of testing to see whether or 
#        not to do a log2 transformation of the data. Since our data is already
#        transformed by VST, we probably don't want the log2 transformation. 
#        When I run the code, "LogC" is FALSE, so I don't think that the log2
#        transformation is actually done. Should probably delete this block?
# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
   ex[which(ex <= 0)] <- NaN
   exprs(gset) <- log2(ex) 
}

# Sean: The design matrix here doesn't take pairing into account I don't think.
#        Maybe we can try running our code to see what the difference is with
#        and without pairing. 
# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2,  adjust="BY", sort.by="B", number=Inf)

# How many probes meet the cutoff?
sum(tT$adj.P.Val < 0.05 & abs(tT$logFC) >= 0.5)
