
setwd("J:/MacLabUsers/Claire/Projects/HVE-microarray/Kavita code and data")

########## KAVITA'S CODE, MY NOTES ##############################
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

load("found_genes.Rda") #gene symbols


#create a mapping from hgnc symbol to go id, strand and entrezid
#set up mart
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get entrez id, go id, strand 
annotations<-getBM(attributes=c("hgnc_symbol","entrezgene",
                                "go_id","strand","namespace_1003"),
                   mart = ensembl)

#filter for just the biological process go ids

annotations_bp<-annotations[annotations$namespace_1003=="biological_process",]

#set up data for the allGenes argument in topGO:
# make a vector of 0's and 1's named with the names of genes
#in the gene univserse. 0 means that
#genes is in the universe but not in your found genes.
#1 means it is in both.

#make a vector of the unique symbols for all_genes to use
#for names in the found vector
CLall_genes<-unique(annotations_bp$hgnc_symbol)

#0's and 1's vector:
#which symbols in all_genes overlap with those in found_genes? (logical)
found2 = factor(as.integer(CLall_genes %in% found_genes))


#assign the names from all_genes to the 1's and 0's in found
names(found2) = CLall_genes

#set up data for GO2Genes argument in topGO

#make a list mapping GO ids and their corresponding symbols
GOtoGene<- split(annotations_bp$hgnc_symbol,annotations_bp$go_id)

# (This also makes a list, but lapply, unique didn't work)
#GotoGene<-unstack(annotations, hgnc_symbol~go_id)

#remove duplicates
GOtoGene<- lapply(GOtoGene, unique)  

#Go2Genes argument needs to be a mapping list
topGOdata <- new("topGOdata",description = "Simple session",
                     ontology = "BP",
                     allGenes = found,
                     nodeSize = 5,
                     annot = annFUN.GO2genes,
                     GO2genes = GOtoGene)

