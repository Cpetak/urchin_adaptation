# Urchin local adaptation
This repo was made to cleanly demonstrate how I got from raw NGS data to different sets of loci putatively under positive selection using Fst outlier, Bayenv and reduced nucleotide diversity measures.

## Step 1: From reads to bams
<details>
  <summary>Click to view detailed code</summary>

folders:  
raw sequence data: WGS/all_fastqs, WGS/all_fastqs_backup, WGS/zipped_backups, WGS/sequencing_meta_data  
scripts and temporary files for alignment and quality check: WGS/scripts_files_for_alignment  
output of alignment: WGS/BWA_out  
	
### Checking quality of sequencing data

```
pip install multiqc
spack load fastqc@0.11.7

-----------
#!/bin/bash

yourfilenames=`ls /users/c/p/cpetak/WGS/all_fastqs/18170X*.fastq`

for file in $yourfilenames

do
	fastqc $file -o /users/c/p/cpetak/WGS/fastqc_output/
done
-----------

cd /users/c/p/cpetak/WGS/fastqc_output
multiqc .
```

#### Results:
[Multiqc Report nicely displayed is available here](https://htmlpreview.github.io/?https://github.com/Cpetak/urchin_adaptation/blob/main/images/multiqc_report.html) 
	
OR 
	
[Multiqc Report raw file for download is available here](images/multiqc_report.html)	


### Mapping to the reference genome

```
spack load bwa@0.7.17
spack load samtools@1.10
bwa index GCF_000002235.5_Spur_5.0_genomic.fna

-----------
while read line ; do
        F1=$(cut -d ' ' -f1 <<< $line)
        F2=$(cut -d ' ' -f2 <<< $line)
        echo "$F1 -- $F2"
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "spack load samtools@1.10" >> $FILE
        echo "spack load bwa@0.7.17" >> $FILE
        ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"
        out_name=$(cut -d '.' -f1 <<< $F1)
        echo "bwa mem -t 1 -M $ref /users/c/p/cpetak/WGS/all_fastqs/$F1 /users/c/p/cpetak/WGS/all_fastqs/$F2 | samtools view -S -b > /users/c/p/cpetak/WGS/BWA_out/$out_name.bam" >> $FILE
          sbatch $FILE
          sleep 0.5
          rm $FILE
done < $1
-----------
```
#### Checking mapping statistics
```
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools sort /users/c/p/cpetak/WGS/BWA_out/$line -o /users/c/p/cpetak/WGS/BWA_out/$out_name.sorted.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools rmdup /users/c/p/cpetak/WGS/BWA_out/$line /users/c/p/cpetak/WGS/BWA_out/$out_name.rmdup.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools flagstat /users/c/p/cpetak/WGS/BWA_out/$line | awk 'NR>=6&&NR<=13 {print \$1}' | column -x >> /users/c/p/cpetak/WGS/$out_name.flagstats.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools depth /users/c/p/cpetak/WGS/BWA_out/$line | awk '{sum+=\$3} END {print sum/NR}' >> /users/c/p/cpetak/WGS/$out_name.coverage.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
```

#### Results:
[File with all mapping stats](all_mapping_stats.csv)
In all 3 images below, x axis is the 140 individuals

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/coverage_fig.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/flagstat_fig.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/mapping_stat.png" width="400" />
	 
</details>

## Step 2: from bams to genotype likelihoods and PCA
I used ANGSD to get genotype likelihoods which then I used to create a PCA and look for population structure and possible outliers
<details>
  <summary>Click to view detailed code</summary>
	
folders:  
scripts and files for PCANGSD: WGS/my_pcangsd/my_PCA (for creating PCA of all individuals based on pops), WGS/my_pcangsd/PCA_bias (for quality check PCA)
	
Run this code on all individuals from all populations together for PCA

```
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
```
python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/pcangsd_covmatrix -threads 16
```
```
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
#### Results:
<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_1.png" width="400" />
	
3 individuals seem to be very different from the other 137 individuals. Thus, these 3 were dropped from further analysis. New PCA with 137 individuals (ANGSD was rerun with only 137 individuals):

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_2.png" width="400" />

No clustering by population can be seen.
	
</details>
	
Then, I proceed to use the PCA for quality check.

<details>
  <summary>Click to view results</summary>
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
</details>	
	
## Step 3: Getting per-site Fst values - global

First, run ANGSD to get bcf file with per site information for each individual from all populations.

```bash
./angsd -b /users/c/p/cpetak/WGS/make_vcf/list_for_angsd.txt \ #all 140 individuals
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/make_vcf/all_pop_angsd \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \ # 85% of 140 individuals
-setMinDepthInd 3 \
-skipTriallelic 1 \
-dobcf 1 \
-GL 1 \
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
-SNP_pval 1e-6
```

-> output of interest: all_pop_angsd.bcf

```bash
# convert to vcf
spack load bcftools@1.10.2
bcftools view all_pop_angsd.bcf > all_pop_angsd.vcf
# then to add GT info
vcfglxgt all_pop_angsd.vcf > fixed_all_pop_angsd.vcf
# then to account for missing information
sed -i 's/0\/0\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/NA\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/g' fixed_all_pop_angsd_copy.vcf
# then keep only GT information
awk -v FS="\t" -v OFS="\t" '{for(i=9;i<=NF;i++) {split($i, gt, ":"); $i=gt[1]} print}' fixed_all_pop_angsd_copy.vcf > fixed_all_pop_angsd_copy_onlyGT.vcf
# getting rid of vcf header
cat fixed_all_pop_angsd_copy_onlyGT.vcf | grep -v "#" > vcf_tail.vcf 
```



For LFMM:

filtered using the bayenv 0025 list:

```bash
# getting first column of vcf as chromosome.position
awk -F "\t" '{print $1$2, $0}' vcf_tail.vcf > vcf_tail_idd2
# keeping only lines in vcf that are also in list of filtered positions
awk 'FNR==NR{a[$0];next}($1 in a)' filt_posi_fixed vcf_tail_idd2 > filtered_vcf
# where filt_posi_fixed is the list of loci after filtering, 1 column, chromosome.position format
# results in a file that has 991,430 positions (994,220 only MAF filter), because remember here we also filtered for SNP_pval
```

then converted into LFMM format after removing CAP and FOG outliers:

```bash
# take out first 10 columns - info about position, and columns 31, 57, 58 which are the outlier individuals
cut --complement -d$'\t' -f1,2,3,4,5,6,7,8,9,30,56,57 filtered_vcf > cut_filtered_vcf
# since file is very big I transposed it using python -> num rows=137, num cols=991,430
with open("cut_filtered_vcf") as f:
    raw = f.readlines()
    data = [el.strip().replace("\n","").split("\t") for el in raw]

def get_col(i,data):
    return [el[i] for el in data]

with open("cut_filtered_vcf_t","w") as f:
    for i in range(len(data[0])):
        f.write(" ".join(get_col(i,data))+"\n")
# make into LFMM format
sed 's/0\/0/0/g' cut_filtered_vcf_t | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' | sed 's/NA/9/g' > lfmm_input
```

Running LFMM:

Getting K

spack load r@3.6.3

```R
#devtools::install_github("bcm-uga/lfmm")
library(lfmm)
library(tseries)

mydata<-read.matrix("lfmm_input", header = FALSE, sep = " ", skip = 0)

#test<-c(0,0,0,2,0,2,2)
#test_pos<-rep(test, each=20)
#mydata2 <- cbind(mydata, test_pos)

# get genotype data, columns are positions, rows are individuals
#Y <- example.data$genotype
Y <- mydata

#principal component analysis (PCA) can reveal some ‘structure’ in the genotypic data. 
#We perfom PCA by using the prcomp function as follows.
pc <- prcomp(Y)
pdf(file = "PCA.pdf")
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
dev.off()
#For the A. thaliana samples, the screeplot indicates that there are around K=6 main components in the data. 
#We will use K=6 latent factors in subsequent analyses of the A. thaliana example.
```

Getting p-values

```R
library(lfmm)
library(tseries)
Y<-read.matrix("lfmm_input", header = FALSE, sep = " ", skip = 0)
#order of vcf: BOD,CAP,FOG,KIB,LOM,SAN,TER
BOD <- rep(0.05026657, each=20)
CAP <- rep(0.07188418, each=19)
FOG <- rep(0.1780976, each=18)
rest <- rep(c(0.003375068,0.00729927,0,0), each=20)
X <- c(BOD, CAP, FOG, rest)

#The R package contains two main functions for estimating the parameters of LFMMs: 
# ridge_lfmm and lasso_lfmm. 
# ridge minimizes sum of squared residuals + lambda * slope squared -> less sensitive to changes in x, lambda is a scaler, the larger lambda the less sensitive y is to x (i.e. flatter the fitted line)
# overall shrinks parameters, making our predictions less sensitive to them
# vs lasso regression: sum of squared residuals + lambda * absolute_value(slope), again making prediction (y) less sensitive to x
# lasso can shrink parameters to 0, ridge can't (when increasing lambda) -> lasso can exclude useless variables from our model
# ridge is more useful when most parameters are useful -> are most of my parameters useful? I think so!
# so I need ridge... let's do both tho and compare results?
# ridge is also the one used in most papers citing the lfmm paper
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 7) #determined from PCA, TODO

#The ridge_lfmm function returns an object that contains the latent variable score matrix U, 
# the latent variable loading matrix U, and B the effect sizes for all SNPs.

## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif") #gif is default, and what others use

pvalues <- pv$calibrated.pvalue 
write.csv(pvalues,"LFMM_ridge_pvalues.csv")
print("created csv with p-vals")

#The histogram of test significance values is expected to be flat, with a peak near zero. 
#A QQ-plot is displayed as follows.

png(file = "QQplot.png")
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
dev.off()

## Manhattan plot
png(file = "Manhattan.png")
plot(-log10(pvalues), 
     pch = 19, 
     cex = .9, 
     xlab = "SNP", ylab = "-Log P",
     col = "black")
dev.off()
```

For FCT I started from filtered_vcf which I copied into a folder named FCT

```bash
#cleaning for outflank
sed -i 's/NA/NA\/NA/g' filtered_vcf
cut --complement -d$' ' -f1 filtered_vcf > cleaned_filtered_vcf
#taking out outliers
cut --complement -d$'\t' -f30,56,57 cleaned_filtered_vcf > cut_filtered_vcf
# put back vcf header
cat vcf_head.vcf cut_filtered_vcf > topped_vcf.vcf
# take out names of outliers in vcf header
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/CAP\_18170X101\_200925\_A00421\_0244\_AHKML5DSXY\_S121\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/FOG\_18170X127\_200925\_A00421\_0244\_AHKML5DSXY\_S147\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/FOG\_18170X128\_200925\_A00421\_0244\_AHKML5DSXY\_S148\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
# removing accidental double tabs
sed 's:\t\t*:\t:g' topped_vcf.vcf > test.vcf
```

Outflank

```R
library(OutFLANK)
library(vcfR)

vcf <- read.vcfR("test.vcf", verbose=FALSE)
ind <- read.table("Pop.txt", header=TRUE)

convertVCFtoCount3 <- function(string){
    a <- as.numeric(unlist(strsplit(string, split = c("[|///]"))))
    odd = seq(1, length(a), by=2)
    a[odd] + a[odd+1]
}
all.vcf.gen <- vcf@gt[,-1]
system.time(gen_table <- matrix(convertVCFtoCount3(all.vcf.gen), ncol=ncol(all.vcf.gen)))

locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")
SNPdata <- t(gen_table)
SNPdata[is.na(SNPdata)] <- 9
k <- max(ind$pop)

FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind$pop)

write.csv(FstDataFrame, file = "results.csv")
```

Where Pop.txt is just a list of 57 1s (first "pop" BOD + CAP + FOG individuals) and 80 2s (second "pop").



## Step 4: Bayenv factor, #TODO

Describe how to get file with pos info for all sites where MAF > 0.025. (bayenv_withpos_0025filter.csv)

## Step 5: Getting per-site Fst values - pair-wise

<details>
  <summary>Click to view detailed code</summary>
</details>

First, I run angsd on each population separately. E.g.

```bash
./angsd -b /users/c/p/cpetak/WGS/BOD_rmdups_jo.txt 
-ref /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-anc /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-out /users/c/p/cpetak/WGS/angsd_new/BOD_angsd_allsites 
-nThreads 16 
-remove_bads 1 
-C 50 
-baq 1 
-minMapQ 30 
-minQ 20 
-minInd 17 # 85% of 20 individuals
-setMinDepthInd 3 
-skipTriallelic 1 
-GL 1 
-doCounts 1 
-doMajorMinor 1 
-doMaf 1 
-doSaf 1 
-doHWE 1
```

Since I have 7 populations, I have 21 possible pairs of populations. For each of the possible pairs:

Note to self!! create new folder, copy all .saf.idx (and other related files) of BOD, TER, SAN, KIB, LOM from angsd_new and CAP, FOG from angsd_noout and do all of the following there

```bash
#!/bin/sh
dir=/users/c/p/cpetak/WGS/angsd_new
while read line ; do #give this script a list of pop pairs displayed as pop1.pop2
    pop1=$(cut -d '.' -f1 <<< $line)
    pop2=$(cut -d '.' -f2 <<< $line)
    echo $pop1
    echo $pop2
    FILE=$(mktemp)
    cat header.txt >> $FILE
    echo "cd /users/c/p/cpetak/WGS/angsd/misc" >> $FILE
    echo "./realSFS ${dir}/${pop1}_angsd_allsites.saf.idx ${dir}/${pop2}_angsd_allsites.saf.idx -P 16 -fold 1 > ${dir}/pairwise_fst/${pop1}_${pop2}_allsites.sfs" >> $FILE #folded option!
    sbatch $FILE
    sleep 0.5
    rm $FILE
done < $1
```

Then for each pair (using TER.BOD as an example)

```bash
cd /users/c/p/cpetak/WGS/angsd/misc
dir=/users/c/p/cpetak/WGS/angsd_new
./realSFS fst index ${dir}/TER_angsd_allsites.saf.idx ${dir}/BOD_angsd_allsites.saf.idx -sfs ${dir}/pairwise_fst/TER_BOD_allsites.sfs -fold 1 -fstout ${dir}/pairwise_fst/TER_BOD_allsites -whichFst 1
```

Finally

```bash
./realSFS fst print /users/c/p/cpetak/WGS/angsd_new/pairwise_fst/TER_BOD_allsites.fst.idx > /users/c/p/cpetak/WGS/angsd_new/pairwise_fst/TER_BOD_allsites.fst
```

To filter by MAF I did the following: keep only rows of the .fst files that are also present in the csv I used for bayenv. (bayenv_onlypos_0025filter.csv)

fixing bayenv and fst files to have correct format:

```bash
sed s/"\.1"/"\.1,"/g bayenv_onlypos_0025filter.csv > bayenv_onlypos_temp.csv
awk -F "," '{print $2"_"$1, $0}' bayenv_onlypos_temp.csv > bayenv_onlypos_temp2.csv
cut -d' ' -f1 bayenv_onlypos_temp2.csv > bayenv_onlypos_temp3.csv
# now it is one column, with pos_chr format, rm intermediate temp files and rename to _cleaned
# repeat with each .fst file too (not exactly as they have different starting format) to match format:
awk '{print $2"_"$1, $0}' CAP_BOD_allsites.fst # example file
# then sort each of the fst files along with the bayenv file
sort -k1 
```

now that the files are all fixed and sorted, we can join each of the fst files with the bayenv file

```bash
# for each fst file
join -1 1 -2 1 sorted_bayenv_onlypos_cleaned.csv ${pop1}_${pop2}_allsites_cleaned.fst > joined_${pop1}_${pop2}_allsites
# IMPORTANT! Make sure that the column you are joining on (in this case the first column) has the same format in both files you are joining! E.g. pos_chr. I chose this IDing instead of chr_pos to avoid sorting issues.
```

Then, created 2 files - one with the list of pop pairs within pH groups, and one with a list of pop pairs between pH groups. For each of these categories I combined the corresponding joined files. 

```bash
cat $(cat between_pairs) > combined_between
cat $(cat within_pairs) > combined_within
```

 I then averaged Fsts coming from within and between pop pairs separately to following way:

```python
import pandas as pd
import matplotlib.pylab as plt
import numpy as np

col_names1=["ID","chr","pos","A", "B"]

df = pd.read_csv('combined_between',names=col_names1, sep="\s")
df["pos"]=df.pos.astype('int64', copy=False)

df["fst"]=df["A"]/df["B"]
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df['fst'] = df['fst'].fillna(0)

df = df.drop(["ID","A","B"], 1)

ndf=df.groupby(["chr","pos"]).mean()

ndf.to_csv('average_between.csv') # do same with within
```

Sometimes the average Fst is slightly negative, let's round that to 0

```bash
awk -F, '$3+0<0{$3=0}1' average_between.csv | sed 's/ /,/g' > average_between_cropped.csv
```

Next, we need to subtract average_within from average_between. I used the following code:

```python
import pandas as pd
import numpy as np

col_names1=["chr","pos", "fst"]

betw = pd.read_csv('average_between_cropped.csv', sep=",")
withi = pd.read_csv('average_within_cropped.csv', sep=",")

betw["pos"]=betw.pos.astype('int64', copy=False)
withi["pos"]=withi.pos.astype('int64', copy=False)
betw["fst"] = pd.to_numeric(betw["fst"], downcast="float")
withi["fst"] = pd.to_numeric(withi["fst"], downcast="float")

cdf=betw.merge(withi, on=["chr","pos"], how="outer")
cdf = cdf.fillna(0)

cdf["final_fst"]=cdf["fst_x"]-cdf["fst_y"]

cdf.to_csv("subbed.csv")
```

subbed.csv contains the final fst value for each pos.

## Step 6: Getting pair-wise global Fst values

<details>
  <summary>Click to view detailed code</summary>
</details>

Method 1:

```bash
./realSFS fst stats ${dir2}/pairwise_fst_cleaned/${pop1}_${pop2}_allsites.fst.idx > ${dir2}/pairwise_fst_cleaned/${pop1}_${pop2}_global.fst
```

-> output: unweighted and weighted (second preferred) global Fst for each pair. However, these are global Fst values prior to filtering based on MAF. So,

Method 2:

```python
import pandas as pd
import matplotlib.pylab as plt
import numpy as np

col_names1=["ID","chr","pos","A", "B"]

df = pd.read_csv('joined_CAP_BOD_allsites',names=col_names1, sep="\s")
df["pos"]=df.pos.astype('int64', copy=False)
#df["A"]=df.A.astype('int64', copy=False)
#df["B"]=df.B.astype('int64', copy=False)

df["fst"]=df["A"]/df["B"]
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df['fst'] = df['fst'].fillna(0)

df = df.drop(["ID","A","B"], 1)

print(df["fst"].mean())
```

Results are in pairwise_global_fst.csv file. My method on the filtered dataset and the realsfs on the original dataset generated pairwise global Fst values that significantly positively correlated (p<0.001, adjusted r-squared = 0.526), and pairwise global Fst values did not correlate with distance between the populations, no matter which method was used, thus there was no isolation by distance in this WGS dataset.

Realsfs vs distance: adj r squared: 0.024, p-val = 0.238

Averaged vs distance: adj r squared: -0.053, p-val = 0.969

Insert 3 images here!

Generated by: https://colab.research.google.com/drive/10XBxj0d1lb9Wp0bcg5VX3XcXvBMolz34?usp=sharing

-> compared to bayenv cov matrix coming!

-> bootstrapped and got outliers!

Bootstrap 99th percentile cut off: 0.02607645058649987

awk -F "," ' $6 >= 0.02607645058649987 ' subbed.csv > pair_fst_outs

 awk -F "," ' $10 >= 0.23131785559999998 ' reps_combined_bayenv.csv > bayenv_outs
