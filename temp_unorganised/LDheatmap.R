install.packages("LDheatmap")

library(LDheatmap)

mat <- scan('resultLD2.ld') # matrix of pairwise LDs, output of PLINK
mat <- matrix(mat, ncol = 100, byrow = TRUE)
locs <- scan('locs_info') # vector of pos info for my SNPs

LDheatmap(mat, genetic.distances = locs, color=heat.colors(20))

