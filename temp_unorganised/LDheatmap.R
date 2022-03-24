install.packages("LDheatmap")

library(LDheatmap)

mat <- scan('resultLD_heat.ld') # matrix of pairwise LDs, output of PLINK
mat <- matrix(mat, ncol = 499, byrow = TRUE)
locs <- scan('locs_info') # vector of pos info for my SNPs

LDheatmap(mat[1:180,1:180], genetic.distances = locs, color=heat.colors(20), add.map=FALSE)
