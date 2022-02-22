install.packages("BiocManager")
BiocManager::install()
BiocManager::install("topGO")
BiocManager::install("ALL")
BiocManager::install("hgu95av2.db")

library(topGO)
library(ALL)

# The geneList data is based on a differential expression analysis 
# of the ALL(Acute Lymphoblastic Leukemia) dataset that was 
# extensively studied in the literature on microarray analysis 
# Chiaretti, S., et al. (2004). Our toy example contains just a 
# small amount, 323, of genes and their corresponding p-values.

data(ALL)
data(geneList)

# The next data one needs are the gene groups itself, the GO terms 
# in our case, and the mapping that associate each gene with one or 
# more GO term(s). The information on where to find the GO annotations 
# is stored in the ALL object and it is easily accessible.

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

# The function assumes that the provided argument is a named vector of p-values. 
# With the help of this function we can see that there are 50 genes with a raw 
# p-value less than 0.01 out of a total of 323 genes.

sum(topDiffGenes(geneList))

# We now have all data necessary to build an object of type topGOdata. This object 
# will contain all gene identifiers and their scores, the GO annotations, the GO 
# hierarchical structure and all other information needed to perform the desired 
# enrichment analysis.

sampleGOdata <- new("topGOdata", description = "Simple session",
                    ontology = "BP", #character string specifying the ontology of interest (BP, MF or CC)
                    allGenes = geneList, #named vector of type numeric or factor. The names attribute contains the genes identifiers. The genes listed in this object define the gene universe.
                    geneSel = topDiffGenes, # function to specify which genes are interesting based on the gene scores
                    nodeSize = 10, # This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated genes (after the true path rule is applied).
                    annot = annFUN.db, 
                    affyLib = affyLib) # function which maps genes identifiers to GO terms.

#  two types of test statistics: Fisher’s exact test which is based on gene counts, 
# and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores. 
# classical enrichment analysis by testing the over-representation of GO terms within
# the group of differentially expressed genes.

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

# other results visualisation in documentation.
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, topNodes = 10)

######

# example of using custom GO map

geneID2GO <- readMappings("geneid2go.map") # Map

geneNames <- names(geneID2GO) # All gene names

myInterestingGenes <- sample(geneNames, length(geneNames) / 10) # list of 10 random geneNames
geneList <- factor(as.integer(geneNames %in% myInterestingGenes)) # mask of 0 and 1 if geneName is interesting
names(geneList) <- geneNames # geneList but annotated with the gene names
str(geneList)

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
allRes

# my stuff
# MIGHT REPEAT WITH OTHER GO MAPPING, FROM ECHINOBASE

geneID2GO <- readMappings("GO_mapping_topGO")
geneNames <- names(geneID2GO)

myInterestingGenes <- read.csv("uniprotIDs_fst_2pops.txt", header = FALSE)
#myInterestingGenes <- as.list(myInterestingGenes["V1"])
intgenes <- myInterestingGenes[, "V1"]
geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
names(geneList) <- geneNames # geneList but annotated with the gene names

GOdata <- new("topGOdata", 
              ontology = "CC", # ontology of interest (BP, MF or CC)
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") # these are the options I'll be using! checked!
resultFisher

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10) # top 10 enriched terms
allRes

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all') # R version error :(

# with echinobase info
# fucking meaningless since echinobase only had 351 genes associated with GO terms
geneID2GO <- readMappings("GO_mapping_topGO_echinoGENEPAGE")
geneNames <- names(geneID2GO)

myInterestingGenes <- read.csv("GENEPAGE_fst_2pops.txt", header = FALSE)
#myInterestingGenes <- as.list(myInterestingGenes["V1"])
intgenes <- myInterestingGenes[, "V1"]
geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
names(geneList) <- geneNames # geneList but annotated with the gene names

GOdata <- new("topGOdata", 
              ontology = "BP", # ontology of interest (BP, MF or CC)
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") # these are the options I'll be using! checked!
resultFisher

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10) # top 10 enriched terms
allRes
