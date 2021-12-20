remotes::install_github("JimWhiting91/genotype_plot")
remotes::install_github("spflanagan/gwscaR")

library(GenotypePlot)
library(gwscaR) # for merging vcfs
library(vcfR)

vcf1 <-parse.vcf("BOD_18170X61_200925_A00421_0244_AHKML5DSXY_S81_L002_R1_001.rmdup.bam_ecm3.vcf")
vcf2 <-parse.vcf("BOD_18170X62_200925_A00421_0244_AHKML5DSXY_S82_L002_R1_001.rmdup.bam_ecm3.vcf")

comb_vcf <- combine.vcfs(vcf1,vcf2, vcf.name = "merge.vcf") #ignores pos only 1 indv is polymorphic

new_plot <- genotype_plot(vcf= system.file("example.vcf.gz", package = "GenotypePlot"),   
                          chr    = 1,                                       # chr or scaffold ID
                          start  = 11700000,                                # start of region
                          end    = 11800000,                                # end = end of region
                          popmap = our_popmap,                              # population membership
                          cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
                          snp_label_size = 10000,                          # breaks for position labels, eg. plot a position every 100,000 bp
                          colour_scheme=c("#FCD225","#C92D59","#300060"),   # character vector of colour values
                          invariant_filter = TRUE)                         # Filter any invariant sites before plotting