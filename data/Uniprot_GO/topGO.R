library(topGO)
setwd("~/Desktop/urchin_adaptation/data/Uniprot_GO")

geneID2GO <- readMappings("GO_mapping_topGO") # uniprot to GO mapping
geneNames <- names(geneID2GO)

myInterestingGenes <- read.csv("temp_files_for_GO/uniprotIDs_prom_locs.txt", header = FALSE) # list of interesting genes, output of LOC to uniprot mapping
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

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 3) # top 10 enriched terms
allRes
write.csv(allRes,"prom_locs_results_CC.csv", row.names = FALSE)


#showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all') 
