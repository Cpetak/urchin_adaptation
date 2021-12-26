# Step 5: running PCAngsd

folders:

WGS/my_pcangsd/persite

Run similar Angsd code for all individuals as before, but less strict (-setMinDepthInd 3 instead of 4) as follows:

```bash
./angsd -b /users/c/p/cpetak/WGS/all_rmdups_noout.txt # all individuals except for the 3 outliers
-ref /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-anc /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-out /users/c/p/cpetak/WGS/angsd_noout/allpopsdepth3_angsd_polysites 
-nThreads 16 
-remove_bads 1 
-C 50 
-baq 1 
-minMapQ 30 
-minQ 20 
-minInd 119 
-setMinDepthInd 3 
-skipTriallelic 1 
-GL 1 
-doCounts 1 
-doMajorMinor 1 
-doMaf 1 
-doGlf 2 
-SNP_pval 1e-6
```

Followed by pcangsd as before:

```bash
python -u /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/my_pcangsd/persite/allpopsdepth3_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/my_pcangsd/persite/PCangsd_selection -selection -sites_save -threads 64
```

TODO repeat with -pcadapt instead of -selection

TODO change minMAF? Default 0.05



Tutorial I was following: http://www.popgen.dk/software/index.php/PCAngsdTutorial

```R
library(RcppCNPy) # Numpy library for R

## function for QQplot
qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
  legend("topleft",paste("lambda=",lambda))
  }

### read in seleciton statistics (chi2 distributed)
s<-npyLoad("PCangsd_selection.selection.npy")

## make QQ plot to QC the test statistics
pdf(file="qqplot.pdf")
qqchi(s)
dev.off()

# convert test statistic to p-value
pval<-1-pchisq(s,1)

## read positions (hg38)
p<-read.table("PCangsd_selection.sites",colC=c("factor","factor","integer"),sep="_")

names(p)<-c("nothing","chr","pos")

write.csv(pval,"pcangsd_results.csv")

## make manhatten plot
pdf(file="manhattenplot.pdf")
plot(-log10(pval),col=p$chr,xlab="Chromosomes",main="Manhattan plot")
dev.off()
```

"The Bonferroni correction compensates for that increase by testing each individual hypothesis at a significance level of α / m where α is the desired overall alpha level and m is the number of hypotheses"

I have 2,119,846 SNPs, alpha = 0.05, -> 2.35x10-8

664 outliers

### Results

Visualised PCA again, code in [here](code/pca.R)

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/pcangsd_madeforpersite.png" width="400" />

No clustering by pop, no clustering my North, Middle, South, no clustering by pH, none of the first 5 PC combinations.