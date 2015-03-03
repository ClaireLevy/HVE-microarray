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
# this is monstrous, but it works. Note that the order is important: because 
# "." gets interpreted as meaning "any character", D1.500 matches both D1.500
# and D14.500. So, D14.500 has to come first so that it will be matched and not
# reach the test for D1.500
# can be verified by comparing levels(factor(df$Condition)) with rename(that)
rename <- function(data) {
   ifelse(str_detect(data, "D14.500"), "D14.500", 
      ifelse(str_detect(data, "D4.500"), "D4.500", 
         ifelse(str_detect(data, "D7.500"), "D7.500", 
            ifelse(str_detect(data, "D1.500"), "D1.500", 
               ifelse(str_detect(data, "D14.50"), "D14.50", 
                  ifelse(str_detect(data, "D4.50"), "D4.50", 
                     ifelse(str_detect(data, "D7.50"), "D7.50", 
                        ifelse(str_detect(data, "D1.50"), "D1.50", 
                           "x"
                        )
                     )
                  )
               )
            )
         )
      )
   )
}

log2FC <- mutate(log2FC, Condition = rename(Condition))
rawFC <- mutate(rawFC, Condition = rename(Condition))
FDR <- mutate(FDR, Condition = rename(Condition))
upDown <- mutate(upDown, Condition = rename(Condition))

################################################################################
# merge them
################################################################################
longForm <- merge(log2FC, rawFC)
b <- merge(FDR, upDown)
longForm <- merge(longForm, b)
rm(FDR, b, log2FC, raw, rawFC, upDown, rename)