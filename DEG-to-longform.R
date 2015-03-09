require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
raw <- read.csv("J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\differentiallyExpressedGenes\\ex_vivo_HVE_Tenofovir_Summary_results.csv")

# strategy: split dataframe into four data frames, melt them, merge 

################################################################################
# split into four
################################################################################

log2FC <- select(raw, Probe.ID, TargetID, ENTREZ_GENE_ID, contains("log2"))
rawFC <- select(raw, Probe.ID, TargetID, ENTREZ_GENE_ID, contains("Symmetrical"))
FDR <- select(raw, Probe.ID, TargetID, ENTREZ_GENE_ID, contains("FDR"))

# some of the columns containing UP DOWN or FALSE are duplicated, e.g., 
# identical(raw$D14.500vNT.DEG, raw$D14.500vNT.DEG.1)
upDown <- select(raw, Probe.ID, TargetID, ENTREZ_GENE_ID, contains("DEG"))
upDown <- select(upDown, -contains("DEG.1"), -contains("ANY"))

################################################################################
# melt them
################################################################################
log2FC <- melt(log2FC, id.vars = c("Probe.ID", "TargetID", "ENTREZ_GENE_ID"), 
               variable.name = "Condition", value.name = "Log2_fold_change")

rawFC <- melt(rawFC, id.vars = c("Probe.ID", "TargetID", "ENTREZ_GENE_ID"), 
              variable.name = "Condition", value.name = "Symmetrical_raw_fold_change")

FDR <- melt(FDR, id.vars = c("Probe.ID", "TargetID", "ENTREZ_GENE_ID"), 
            variable.name = "Condition", value.name = "FDR_adjusted_p_value")

upDown <- melt(upDown, id.vars = c("Probe.ID", "TargetID", "ENTREZ_GENE_ID"), 
               variable.name = "Condition", value.name = "Direction")

################################################################################
# rename them
################################################################################
# Note that the order is important:  D1 matches both D1 and D14. So, D14 has to 
# come first so that it will be matched and not reach the test for D1
# can be verified by comparing levels(factor(df$Condition)) with rename(that)
days <- function(data) {
  ifelse(str_detect(data, "D14"), 14, 
         ifelse(str_detect(data, "D4"), 4, 
                ifelse(str_detect(data, "D7"), 7, 
                       ifelse(str_detect(data, "D1"), 1, 0
                       )
                )
         )
  )
}

# same order comment
concs <- function(data) {
  ifelse(str_detect(data, "500"), 500, 50)
}

log2FC <- mutate(log2FC, Day = days(Condition), 
                 Concentration = concs(Condition)) %>%
  select (-(Condition))
rawFC <- mutate(rawFC, Day = days(Condition), 
                Concentration = concs(Condition)) %>%
  select (-(Condition))
FDR <- mutate(FDR, Day = days(Condition), 
              Concentration = concs(Condition)) %>%
  select (-(Condition))
upDown <- mutate(upDown, Day = days(Condition), 
                 Concentration = concs(Condition)) %>%
  select (-(Condition))

################################################################################
# merge them
################################################################################
longForm <- merge(log2FC, rawFC)
b <- merge(FDR, upDown)
longForm <- merge(longForm, b)
rm(FDR, b, log2FC, raw, rawFC, upDown, concs, days)