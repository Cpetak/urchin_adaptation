install.packages("BiocManager")
BiocManager::install()
BiocManager::install("topGO")
BiocManager::install("ALL")
BiocManager::install("hgu95av2.db")

library(topGO)
library(ALL)
#library(affy)

data(ALL)
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = topDiffGenes,nodeSize = 10,annot = annFUN.db, affyLib = affyLib)
