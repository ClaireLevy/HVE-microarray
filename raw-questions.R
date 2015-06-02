################# Probe retention through pre-filtering #######################
### Are all of the probes that Lamar found to be differentially expressed 
### retained after pre-filtering?

# run code in "working with raw data.R" to get the following objects:
#     allDays.N.Q = all probes
#     expressed = probe detected in all 3 cell lines for at least one condition
#     dataMatrixSDfilter = expressed AND row SDs >= 0.2

# get vector of 7832 probes from Lamar's Cyber-T analysis
source("H:\\R\\HVE-microarray\\DEG-to-long-form.R") # get longForm 
probes <- unique(longForm$Probe.ID)

# compare probes retained at various stages of prefiltering
sum(probes %in% rownames(exprs(allDays.N.Q))) # 7832
sum(probes %in% rownames(exprs(expressed))) # 7825 (lost 7)
sum(probes %in% rownames(dataMatrixSDfilter)) # 7756 (lost 69 more)

# total loss in prefiltering = 76

######################## Design matrix order ##################################
### Does the order of the design matrix order matter? 

# compare the following two design matrices

# designA
#     CellLine    Condition
#        1           Ctrl
#        2           Ctrl
#        3           Ctrl
#        1           Drug
#        2           Drug
#        3           Drug
targets <- data.frame("CellLine" = c("HVE1", "HVE2", "HVE3", "HVE1", "HVE2",
      "HVE3"),
   "Treatment" = c("control", "control", "control", "drug", "drug", "drug"))
CellLine <- factor(targets$CellLine, levels=c("HVE1", "HVE2", "HVE3"))
Treatment <- factor(targets$Treatment, levels=c("control", "drug"))
designA <- model.matrix(~0 + CellLine + Treatment)

# designB
#     CellLine    Condition
#        1           Ctrl
#        1           Drug
#        2           Ctrl
#        2           Drug
#        3           Ctrl
#        3           Drug
targets <- targets[c(1, 4, 2, 5, 3, 6), ]
CellLine <- factor(targets$CellLine, levels=c("HVE1", "HVE2", "HVE3"))
Treatment <- factor(targets$Treatment, levels=c("control", "drug"))
designB <- model.matrix(~0 + CellLine + Treatment)

library(limma)
load("J:/MacLabUsers/Claire/Projects/HVE-microarray/microarrayData/dataMatrixSDfilter.Rda")# data filtered by SD and detection
dataA <- dataMatrixSDfilter[ , c("HVE_A4", "HVE_A8", "HVE_A12", 
   "HVE_C4", "HVE_C8", "HVE_C12")]
dataB <- dataMatrixSDfilter[ , c("HVE_A4", "HVE_C4", "HVE_A8", 
   "HVE_C8", "HVE_A12", "HVE_C12")]

fitA <- eBayes(lmFit(dataA, design = designA))
ttA <- topTable(fitA, coef = "Treatmentdrug", adjust = "BH", number = 10380)
sum(ttA$adj.P.Val < 0.05)
# -> 608 probes

fitB <- eBayes(lmFit(dataB, design = designB))
ttB <- topTable(fitB, coef = "Treatmentdrug", adjust = "BH", number = 10380)
sum(ttB$adj.P.Val < 0.05)
# -> 608 probes

all(rownames(ttA) == rownames(ttB)) # TRUE
isTRUE(all.equal(ttA, ttB)) # TRUE
# note: isTRUE(all.equal()) is similar to identical, but allows small variation
# in numbers (difference < 1.5e-8) to account for imprecision in how R stores
# numbers

# note2: can't test fit objects for equality directly because the fit$design
# differs between them
lapply(names(fitA), FUN = function(name) 
      isTRUE(all.equal(fitA[[name]], fitB[[name]])))
# -> TRUE for all but fit$design and fit$qr