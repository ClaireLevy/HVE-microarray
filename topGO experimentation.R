library(dplyr)


#make a vector of the ids to annotate

source("DEG-to-long-form.R")
up50Day7<-longForm%>%
  filter(Day=="7", Concentration=="50",Direction=="UP")

library(topGO)
library(GOstats)
library(GO.db)
library(annotate)
library(genefilter)
library(RColorBrewer)
library(xtable)
library(Rgraphviz)
library(illuminaHumanv4.db)
library(biomaRt)


#remove na's and select unique entrez ids only
myGenes<-na.omit(up50Day7)
myGenesEntrez<-myGenes$ENTREZ_GENE_ID
myGenesEntrez<-unique(myGenesEntrez)
###########################################################

#topGO

#get GO ids for all the genes in the array universe
#this is a list of entrez ids and under each are the go numbers
#associated with those entrez ids.

#different options for how to get the GO's
#annFUN.db(whichOnto, feasibleGenes = NULL, affyLib)
#annFUN.org(whichOnto, feasibleGenes = NULL, mapping, ID = "entrez") 
#annFUN(whichOnto, feasibleGenes = NULL, affyLib)

#The functions annFUN.gene2GO and annFUN.GO2genes are used 
#when the user provide his own annotations either as a
#gene-to-GOs mapping, either as a GO-to-genes mapping.

#The annFUN.org function is using the mappings from the 
#"org.XX.XX" annotation packages. The function supports 
#different gene identifiers.

#The annFUN.file function will read the annotations of
#the type gene2GO or GO2genes from a text file.

geneID2GO<-inverseList(annFUN.org("BP", feasibleGenes = NULL,
                      mapping="org.Hs.eg.db", ID = "entrez"))


#geneNames are the entrez ids associated with GO numbers
geneNames<-names(geneID2GO)


#make a list of 1's and 0's to show which genes in the 
#universe overlap with my interesting genes.

geneList<-factor(as.integer(geneNames %in% myGenesEntrez))
#geneList stays the same length as geneNames since the entries
# not actually filtered here, just marked 0 or 1
#so geneNames will have the correct number of names to use for the 
#geneList

#make geneList a named list, using the names from geneNames
names(geneList)<-geneNames


#allGenes is the list showing which of my genes are in the universe
#(1) and which are not(0)

## there are three annotation functions available:
##1. annFUN.db  -- used for bioconductor annotation chips
##2. annFUN.gene2GO  -- used when you have mappings from each gene to GOs
##3. annFUN.GO2genes -- used when you have mappings from each GO to genes

#I have gene 2 GOs in the geneID2GO list that I made so I will
#use option 2.
#gene2GO is the mapping that shows geneIDs and assoc GO numbers
#in the correct list format.


GOdata<-new("topGOdata", ontology = "BP", allGenes = geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO)




head(termStat(GOdata))

#stats for one GO id
GOID<-"GO:0000003"
gene.universe<-genes(GOdata)
go.genes<-genesInTerm(GOdata,GOID)[[1]]
sig.genes<-sigGenes(GOdata)

my.group<-new("classicCount",testStatistic=GOFisherTest,
              name="fisher",allMembers=gene.universe,
              groupMembers=go.genes,
              sigMembers=sig.genes)

contTable(my.group)

runTest(my.group)
test.stat<-new("classicCount",
               testStatistic = GOFisherTest,
               name = "Fisher test")
resultFisher<-getSigGroups(GOdata, test.stat)

#I guess this is for all the GOdata?
resultFis<-runTest(GOdata,
                   algorithm="classic",statistic="fisher")

pvalFis<-score(resultFis)
head(pvalFis)
#you can see the pval for GO:0000003 is the same as what I got
#above when I got the stats for just that GO id.

#plot plot plot

hist(pvalFis,50,xlab="p-values")

#an ordered table showing results from different tests
#I am only showing the classic one since I haven't done the others
#gives top 20 results incl term names.

allRes <- GenTable(GOdata, classic = resultFis,
                   orderBy = "weight",
                   ranksOf = "classic", topNodes = 20)
allRes

#show the top 5 significant results on a GO plot
#based on the results from the fisher test.


showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
