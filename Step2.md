# Step 2: from bams to genotype likelihoods and PCA

I used ANGSD to get genotype likelihoods which then I used to create a PCA and look for population structure and possible outliers

folders:  
scripts and files for PCANGSD: WGS/my_pcangsd/my_PCA 
	
Run this code on all individuals from all populations together for PCA

```bash
cd /users/c/p/cpetak/WGS/angsd

ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna" 
# latest version of the reference genome downloaded from NCBI on the 7th of October 2020
# angsd version: 0.933-102-g7d57642 (htslib: 1.11-9-g2264113) build(Oct 16 2020 18:14:45)

./angsd -b /users/c/p/cpetak/WGS/all_rmdups_jo.txt \
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \ # 85% of all individuals (140)
-setMinDepthInd 4 \ # note that later we'll use 3 here, filtering is stricter for now to reduce data for PCA
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \ # gives us the Beagle format which will be used by pcangsd
-SNP_pval 1e-6
	
```

```bash
python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/pcangsd_covmatrix -threads 16
```

```R
#R
C <- as.matrix(read.table("pcangsd_covmatrix.cov"))
ids <- read.table("~/Downloads/pca_pops.txt") #text file with 20 lines of the single word BOD, then 20 lines of CAP etc in the order they appeared in all_rmdups_jo.txt
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
library(tidyverse)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()
```

### Results:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_1.png" width="400" />
	
3 individuals seem to be very different from the other 137 individuals. Thus, these 3 were dropped from further analysis. New PCA with 137 individuals (ANGSD was rerun with only 137 individuals):

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_2.png" width="400" />

No clustering by population can be seen.
	

## Then, I proceed to use the PCA for quality check.

folders:  
scripts and files for PCANGSD: WGS/my_pcangsd/PCA_bias

I checked if there is any clustering by coverage.

Histogram of average coverage per individual:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/hist_coverage.png" width="400" />
	

```
#R
covdata <- read.table("covs.txt") #this is without outliers, average coverage per individual
covdata$V2 <- ifelse(covdata$V1<6.5, "little", "lot")

ids <-covdata

#ggplot
df <- data.frame(pop = ids$V2, PC1 = e$vectors[,1], PC2 = e$vectors[,2])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()
```

PCA of average coverage:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_cov.png" width="400" />
     
Again, no clustering is visible.
	