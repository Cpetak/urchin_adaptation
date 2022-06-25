library(topGO)
setwd("~/Desktop/urchin_adaptation/data/Uniprot_GO")

geneID2GO <- readMappings("GO_mapping_topGO") # uniprot to GO mapping
geneNames <- names(geneID2GO)

dotopgo <- function(myfile, gokind, outname){
  myfilename<-paste("temp_files_for_GO/", myfile, sep="")
  print(myfilename)
  myInterestingGenes <- read.csv(myfilename, header = FALSE) # list of interesting genes, output of LOC to uniprot mapping
  intgenes <- myInterestingGenes[, "V1"]
  geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
  names(geneList) <- geneNames # geneList but annotated with the gene names
  
  GOdata <- new("topGOdata", 
                ontology = gokind, # ontology of interest (BP, MF or CC)
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO)
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") # these are the options I'll be using! checked!
  
  allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = sum(resultFisher@score < 0.01)) # top 10 enriched terms
  
  write.csv(allRes,outname, row.names = FALSE)
}

for (j in list("all_locs","prom_locs","non_prom_locs")){
  for (i in list("BP","MF","CC")){
    inputname<-paste("uniprotIDs_",j,".txt", sep="")
    dotopgo(inputname, i, paste(j,"_results_",i,"_888_5000.csv", sep=""))
  }
}




#showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all') 
