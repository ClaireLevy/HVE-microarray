require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
source("getAllDavid.R")
source("getAllInnateDB.R")
source('subsetToOverlappingGoTerms.R')

###############################################################################
## Heat map of genes associated with GO term 005615 (extracellular space)
##
## 1. Read in data from DAVID and InnateDB
## 2. Find the days/concentrations/directions where the term is enriched
## 3. Get all of the genes associated with the term in those conditions,
##    combining all genes together. 
## 4. Get the gene expression data for those genes on all days/conditions
## 5. Make the heat map (sorting the genes by the number of days on which they
##    were changed; i.e. genes that changed on all 4 days together, 3 days 
##    together, etc.)

## Step 1: Read in data from DAVID and InnateDB
allD50 <- getAllDavid(50)
allD500 <- getAllDavid(500)
allInnate50 <- getAllInnateDB(50)
allInnate500 <- getAllInnateDB(500)

head(allD50)
# > Source: local data frame [6 x 16]
# > 
# > Category                                           Term Count Percentage
# > 1  PANTHER_MF_ALL             MF00218:Calmodulin related protein     3   20.00000
# > 2   GOTERM_BP_FAT               GO:0008544~epidermis development     3   20.00000
# > 3   GOTERM_BP_FAT                GO:0007398~ectoderm development     3   20.00000
# > 4  PANTHER_MF_ALL         MF00188:Select calcium binding protein     3   20.00000
# > 5 SP_PIR_KEYWORDS                                     Ichthyosis     2   13.33333
# > 6    KEGG_PATHWAY hsa04070:Phosphatidylinositol signaling system     2   13.33333
# > Variables not shown: PValue (dbl), Genes (chr), List.Total (int), Pop.Hits (int),
# > Pop.Total (int), Fold.Enrichment (dbl), Bonferroni (dbl), Benjamini (dbl), FDR (dbl),
# > Day (dbl), direction (chr), Pathway.Id (chr)

## Step 2: Find the days/concentrations/directions where the term is enriched
any(str_detect(allD50$Term, "GO:0005615"))         # yes
any(str_detect(allD500$Term, "GO:0005615"))        # no
any(str_detect(allInnate50$Term, "GO:0005615"))    # yes
any(str_detect(allInnate500$Term, "GO:0005615"))   # yes

selectColumnsOfInterest <- function(df) {
   z <- df %>%
      filter(str_detect(Term, "GO:0005615"), 
         PValue < 0.05)
   
   if("Gene.Symbols" %in% colnames(z)) {
      dplyr::select(z, Day, PValue, Count, direction, Gene.Symbols)
   } else {
      dplyr::select(z, Day, PValue, Count, direction, Genes)
   }
}

selectColumnsOfInterest(allD50) # down d14
selectColumnsOfInterest(allD500) # none
selectColumnsOfInterest(allInnate50) # down d4, 7, 14
selectColumnsOfInterest(allInnate500) # down d4, 7, 14

## Step 3: Get the genes
# convert david genes to symbols
entrez <- str_split(selectColumnsOfInterest(allD50)$Genes, ", ")
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
david50Symbols <- getBM(attributes = c("hgnc_symbol"), filters = "entrezgene",
   values = entrez, mart = mart)[[1]]

innate50Symbols <- str_split(selectColumnsOfInterest(allInnate50)$Gene.Symbols, "; ")
innate500Symbols <- str_split(selectColumnsOfInterest(allInnate500)$Gene.Symbols, "; ")
genes <- unique(unlist(c(david50Symbols, innate50Symbols, innate500Symbols)))

# lots of overlap
length(unique(david50Symbols)) # 52
length(unique(unlist(innate50Symbols))) # 161
length(unique(unlist(innate500Symbols))) # 179
length(unique(unlist(c(innate50Symbols, innate500Symbols)))) # 192
length(genes) # 197

## Step 4: Get the gene expression data for those genes
# get longForm 
source("DEG-to-long-form.R")
source("collapse-probes.R")
z <- longForm %>%
   group_by(Direction, Concentration, Day) %>%
   collapseProbes() %>%
   filter(TargetID %in% genes, Direction %in% c("UP", "DOWN", "FALSE")) 
head(z)
# > Source: local data frame [6 x 8]
# > Groups: Day, Concentration
# > 
# > Day Concentration TargetID ENTREZ_GENE_ID Log2_fold_change
# > 1   1            50    A2ML1         144568      -0.58544610
# > 2   1            50   ABI3BP          25890       0.17863560
# > 3   1            50    ACTA2             59      -0.03381754
# > 4   1            50    ACTG1             71      -0.15087730
# > 5   1            50    ACTN4             81      -0.34227020
# > 6   1            50      ADM            133       0.01285804
# > Variables not shown: Symmetrical_raw_fold_change (dbl), FDR_adjusted_p_value (dbl),
# > Direction (chr)

## Step 5: Make the heatmap
# sort the genes
genesInOrder <- z %>%
   ungroup() %>%
   filter(Concentration == 500) %>%
   # add leading 0 to days < 10 for better ordering (01, 04, 14 not 1, 14, 4)
   mutate(Day = stringr::str_pad(Day, 2, pad = "0")) %>%
   # collapse to one row per target by targetID
   group_by(TargetID) %>%
   summarize(n = sum(Direction != "FALSE"),
      up = sum(Direction == "UP"),
      down = sum(Direction == "DOWN"),
      which = paste(Day[Direction != "FALSE"], collapse = "."), 
      Direction = ifelse(up > 0 & down > 0, "BOTH", 
         ifelse(up > 0, "UP", 
            ifelse(down > 0, "DOWN", "FALSE")))) %>%
   # put in order of the new variables
   arrange(Direction, n, which) %>%
   # extract TargetID column as vector
   .$TargetID

head(genesInOrder)
# > [1] FSTL3  ANXA1  APP    TGFA   IL32   MTHFD2
# > 6699 Levels: 41337 41339 41519 41522 41524 41528 41532 A2ML1 A4GALT AACS ... ZZZ3

toPlot <- z %>%
   mutate(TargetID = factor(TargetID, levels = as.character(genesInOrder)))

# make the plot
p <- ggplot(toPlot, aes(x = factor(Day), y = TargetID))
p + geom_tile(aes(fill = Direction)) + facet_wrap( ~ Concentration) +
   scale_fill_brewer(type = "div", palette = 5) +
   ggtitle("GO:0005615 - extracellular space genes") +
   xlab("Day") + ylab("Gene") + theme(axis.text.y = element_text(size = 4))

ggsave("extracellular-space-genes-up-or-down.pdf", 
   height = 11, width = 8.5, units = "in")



# attempts at using log2 fold change instead of just up or down below
p + geom_tile(aes(fill = Log2_fold_change)) + facet_wrap( ~ Concentration) +
   scale_fill_gradient2(low = "#EF8A62", mid = "#F7F7F7", high = "#67A9CF") +
   ggtitle("GO:0005615 - extracellular space genes") +
   xlab("Day") + ylab("Gene") + theme(axis.text.y = element_text(size = 4))

p + geom_tile(aes(fill = Log2_fold_change)) + facet_wrap( ~ Concentration) +
   scale_fill_distiller(type = "div", palette = 5, space = "rgb") +
   ggtitle("GO:0005615 - extracellular space genes") +
   xlab("Day") + ylab("Gene") + theme(axis.text.y = element_text(size = 4))


# bin continuous data for nicer heatmap
bins <- seq(min(toPlot$Log2_fold_change), -1*min(toPlot$Log2_fold_change), length.out = 11)

toPlot2 <- mutate(toPlot, Discrete_Log2_Fold_Change = 
      bins[findInterval(Log2_fold_change, bins)]) 


p2 <- ggplot(toPlot2, aes(x = factor(Day), y = TargetID))
p2 + geom_tile(aes(fill = factor(Discrete_Log2_Fold_Change))) + facet_wrap( ~ Concentration) +
   scale_fill_brewer(type = "div", palette = 5) + 
   ggtitle("GO:0005615 - extracellular space genes") +
   xlab("Day") + ylab("Gene") + theme(axis.text.y = element_text(size = 4))




# write a file mapping gene symbol to entrez ID for cases of ambiguity
toPlot %>%
   ungroup() %>%
   filter(Day == 1, Concentration == 50) %>%
   mutate(GeneSymbol = TargetID, EntrezGeneId = ENTREZ_GENE_ID) %>%
   # put in reverse order of TargetID to match heatmap
   arrange(desc(TargetID)) %>%
   dplyr::select(GeneSymbol, EntrezGeneId) %>%
   write.csv(
      file = "C:/Users/smhughes/desktop/gene-symbol-to-entrez-id-mapping.csv", 
      row.names = FALSE)