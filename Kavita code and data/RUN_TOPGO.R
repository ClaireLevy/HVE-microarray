
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/Kavita code and data")

########## KAVITA'S CODE, MY NOTES##############################
require(topGO)
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/Kavita code and data")
load("all_genes.Rda") #universe? gene symbols
load("entrez_annot.Rda") #mapping GO to gene symbols
load("found_genes.Rda") #gene symbols
load("found_genes_wang.Rda")#gene symbols


found = factor(as.integer(all_genes %in% found_genes))

#assign the names from all_genes to the 1's and 0's in found
names(found) = all_genes

#why bp only?
PsampleGOdata <- new("topGOdata",description = "Simple session",
                     ontology = "BP",
                     allGenes = found,
                     nodeSize = 5,
                     annot = annFUN.GO2genes,
                     GO2genes = entrez_annot)

PresultClassic = runTest(PsampleGOdata,algorithm = 'classic',
                         statistic =  'fisher')

print(PresultClassic)

#I don't understand what the elim algorithm does, it is used in
 #different context in the vignette
Presultelim = runTest(PsampleGOdata,algorithm = 'elim',
                      statistic =  'fisher')

print(Presultelim)

#showing results for both algorithms
#Rclassic and Relim are column names??
#topNodes is how many to show, here it is rows in PresultClassic score column
#ranksOf shows what the rank of that term was using the Relim algorithm
#orderBy says to order the results by the results of the classic algorithm

Presults.table = GenTable(PsampleGOdata, Rclassic = PresultClassic,
                          Relim = Presultelim,
                          topNodes = length(PresultClassic@score), 
                          ranksOf = "Relim", orderBy = "Rclassic")

#a GenTable is a df so it can be subsetted normally
Presults.table.bh <- Presults.table[as.numeric(Presults.table$Rclassic) < 0.01,]

write.table(Presults.table.bh,"SCRI.GO.tsv",sep = "\t")

#Here is the code for getting genes for a GO term

#genesInTerm 
for (i in 1:nrow(Presults.table)) {
        GOid.of.interest = Presults.table.bh[i,"GO.ID"]
        all.term.genes = genesInTerm(PsampleGOdata,GOid.of.interest)[[1]]
        genes = intersect(found_genes,all.term.genes);
        line = paste(GOid.of.interest,genes,sep = "\t");
        write.table(line,"GenesInGoSCRI.tsv",append = TRUE)
}


######################### MY CODE AND NOTES ################################
require(biomaRt)
require(topGO)
load("all_genes.Rda") #universe? gene symbols
load("found_genes.Rda") #gene symbols
load("found_genes_wang.Rda")#gene symbols


#which in all_genes overlap with those in found_genes (logical)
found = factor(as.integer(all_genes %in% found_genes))

#assign the names from all_genes to the 1's and 0's in found
names(found) = all_genes

#set up mart
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get entrez id, go id, strand for found_genes
annotations<-getBM(attributes=c("hgnc_symbol","entrezgene","go_id","strand"),
             filters = "hgnc_symbol",values= found_genes,
             mart = ensembl)

#make a list mapping GO ids and their corresponding symbols (format for topGO)
GOtoGene<- split(annotations$hgnc_symbol,annotations$go_id)
GOtoGene<- lapply(GOtoGene, unique)  # to remove duplicates

#now  have a new mapping of GO to symbol,can use topGO

topGOdata <- new("topGOdata",description = "Simple session",
                     ontology = "BP",
                     allGenes = found,
                     nodeSize = 5,
                     annot = annFUN.GO2genes,
                     GO2genes = symbolToGO)