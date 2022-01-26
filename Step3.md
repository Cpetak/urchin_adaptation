# Step 3: Getting per-site Fst values - global

folders:  

scripts and files: WGS/make_vcf (this is where I made a filtered large vcf file for all individuals (vcf_head.vcf and vcf_tail.vcf) and subdirectories that contain results of steps that use this vcf)

## Make vcf file

First, run ANGSD to get bcf file with per site information for each individual from all populations.


```bash
./angsd -b /users/c/p/cpetak/WGS/make_vcf/list_for_angsd.txt \ #all 140 individuals!!!
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
head -916 fixed_all_pop_angsd_copy_onlyGT.vcf > vcf_head.vcf
cat fixed_all_pop_angsd_copy_onlyGT.vcf | grep -v "#" > vcf_tail.vcf 
```

Note:  fixed_all_pop_angsd.vcf and fixed_all_pop_angsd_onlyGT.vcf in the make_vcf folder are the files prior the the above 2 steps!!! so ignore them and use vcf_tail.vcf

## Use Outflank for per-site Fst - 7 pops

NOTE: I didn't end up using this version of the code... look below

Folder: 

WGS/make_vcf/using_vcf/my_outflank, small temporary files: /WGS/make_vcf/using_vcf/my_outflank/small_temp

Since the vcf file is huge with many SNPs (unfiltered, all individuals (140), 15,902,843), I divided it up into smaller files, run outflank on each chunk, then concatenated the results as follows:

```bash
split -l 100000 vcf_tail.vcf
# putting back header for each temporary file
#!/bin/bash
pfiles=$(ls | grep "x")
for i in $pfiles
do
        cat vcf_head.vcf $i > topped_${i}.vcf
done
# fixing file format
#!/bin/bash
pfiles=$(ls | grep "topped")
for i in $pfiles
do
   echo $i
   sed -i 's/NA/NA\/NA/g' $i
done
# separate R script for each partial vcf to analyse
#!/bin/sh
pfiles=$(ls | grep "topped")
for i in $pfiles
do
	echo $i
	script_name=${i}script.R
	echo -e "library(OutFLANK)" >> $script_name
	echo -e "library(vcfR)" >> $script_name
	echo -e "vcf <- read.vcfR(\"$i\", verbose=FALSE)" >> $script_name
	echo -e "ind <- read.table(\"Pop.txt\", header=TRUE)" >> $script_name
	echo -e "convertVCFtoCount3 <- function(string){" >> $script_name
	echo -e "\ta <- as.numeric(unlist(strsplit(string, split = c(\"[|///]\"))))" >> $script_name
	echo -e "\todd = seq(1, length(a), by=2)" >> $script_name
	echo -e "\ta[odd] + a[odd+1]" >> $script_name
	echo -e "}" >> $script_name
	echo -e "all.vcf.gen <- vcf@gt[,-1]" >> $script_name
	echo -e "system.time(gen_table <- matrix(convertVCFtoCount3(all.vcf.gen), ncol=ncol(all.vcf.gen)))" >> $script_name
	echo -e "locinames <- paste(vcf@fix[,\"CHROM\"], vcf@fix[,\"POS\"], sep=\"_\")" >> $script_name
	echo -e "SNPdata <- t(gen_table)" >> $script_name
  echo -e "SNPdata[is.na(SNPdata)] <- 9" >> $script_name
	echo -e "k <- max(ind\$pop)" >> $script_name
	echo -e "FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind\$pop)" >> $script_name
	echo -e "write.csv(FstDataFrame, file = \"${i}data.csv\")" >> $script_name
done
#launch each r script as a separate job
while read line ; do
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "Rscript $line" >> $FILE
        sbatch $FILE
        sleep 0.5
        rm $FILE
done < $1
#concat csvs
mv *data.csv csvs
cd csvs
for i in $(ls); do sed '1d' $i > ${i}_fixed; done #fixing csvs
cat *csv_fixed > combined.csv #moved temp files into small_temp
cut -f 2-10 -d , combined.csv | nl -w 1 -p -s , > fixed_combined.csv #reindexing csv
vim fixed_combined.csv -> insert as first line: "","LocusName","He","FST","T1","T2","FSTNoCorr","T1NoCorr","T2NoCorr","meanAlleleFreq" 
#moved to results folder
awk -F, '$3 > 0.1' fixed_combined.csv > fixed_combined_goodhe.csv #getting only sites with He > 0.1
cat fixed_combined_goodhe.csv | cut -d ',' -f2,7 > twocol.csv #keeping only position and FSTNoCorr columns
sort -k 2 -t , -n -r twocol.csv > sorted_twocol.csv #sort by Fst
```

Where Pop.txt is just a list of 1 repeated 20 times, 2, repeated 20 times, etc.

To summarise, there are results of angsd -> vcf of all individuals, no filter, all 15 million sites

In fixed_combined_goodhe.csv: only 2,625,660 as these are He > 0.1

## Filtering vcf_tail.vcf

folder:

WGS/make_vcf/using_vcf/LFMM

filtered vcf_tail.vcf using the bayenv 0025 list* as follows:

```bash
# getting first column of vcf as chromosome.position
awk -F "\t" '{print $1$2, $0}' vcf_tail.vcf > vcf_tail_idd2
# keeping only lines in vcf that are also in list of filtered positions
awk 'FNR==NR{a[$0];next}($1 in a)' filt_posi_fixed vcf_tail_idd2 > filtered_vcf
# where filt_posi_fixed is the list of loci after filtering, 1 column, chromosome.position format
# results in a file that has 991,430 positions (994,220 only MAF filter*), because remember here we also filtered for SNP_pval
```

*This is coming from the filtering steps described in [Code for Filtering steps](https://github.com/Cpetak/urchin_adaptation/blob/main/Filtering_steps.md)

## Use Outflank for per-site Fst - 2 pops

folder:

WGS/make_vcf/using_vcf/FCT

I started from filtered_vcf as described above (MAF filtering), then cleaned for Outflank as follows:

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

```R
library(OutFLANK)
library(vcfR)

FstDataFrame <- read.csv(file = 'results.csv', header=TRUE,row.names=1)

#reduced_df<-FstDataFrame[seq(1,nrow(FstDataFrame),1000),]
reduced_df<-FstDataFrame

pdf("line.pdf")
plot(reduced_df$FST, reduced_df$FSTNoCorr, xlim=c(-0.01,0.3), ylim=c(-0.01,0.3), pch=20)
abline(0,1)
dev.off()

pdf("dots.pdf")
plot(reduced_df$He, reduced_df$FSTNoCorr, pch=20, col="grey")
dev.off()

pdf("hist.pdf")
hist(reduced_df$FSTNoCorr[reduced_df$He>0.1],xlim=c(0,0.3), breaks=50)
dev.off()
```

```R
library(OutFLANK)
library(vcfR)

FstDataFrame <- read.csv(file = 'results.csv', header=TRUE,row.names=1)

print("read file")

k=2
q=0.05

outlier <- OutFLANK(FstDataFrame, NumberOfSamples=k) #investigate options

write.csv(outlier, file = "outliers_from_outflank.csv")

print("created outlier")

pdf("outflank.pdf")
OutFLANKResultsPlotter(outlier, withOutliers = TRUE,NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) #investigate options
dev.off()

print("created outflank.pdf")

pdf("p_hist.pdf")
hist(outlier$results$pvaluesRightTail)
dev.off()

num_out <- sum(outlier$results$qvalues<q, na.rm=TRUE)

if (num_out > 0) {
print("there are outliers:")
print(num_out)
pdf("outliers.pdf")
plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
    points(outlier$results$He[outlier$results$qvalues<q], outlier$results$FST[outlier$results$qvalues<q], pch=21, col="blue")
dev.off()

print("created outliers.pdf")

top_candidates <- outlier$results$qvalues<q & outlier$results$He>0.1
topcan <- outlier$results[top_candidates,]

write.csv(topcan, file = "top_fst.csv")
}

print("all done")
```

### Results

For all figures, see [this folder](https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT). Note: Don't try to open line and dot files saved as pdf, they are huge and will likely freeze your computer :)

I calculated Fsts 4 ways: for 2 and 7 populations, and for sites with random filtering (same number of SNPs in the end but selected at random, regardless of MAF) and filtering as described above (considering MAF, etc), to see if my way of filtering influenced the overall distribution of Fsts

All line and dot plots looked good. A description of what these plots show: http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html

A representative line plot (for 2 pop, normal filtering): 

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/line_2pop.jpg" width="400" />

SNPs with missing a lot of data would have an elevated value of uncorrected Fst relative to corrected Fst. SNPs like this should be removed before running OutFLANK.

A representative dots plot (for 7 pop, normal filtering):

 <img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/dots_7pop.jpg" width="400" />

This is just meant to show how for low He SNPs Fst is inflated. SNPs for which He < 0.1 were removed in outlier calculations (OutFlank function, default). Note how there aren't any positions He < 0.05, as there sites were filtered out already.

Now let's look at the distribution of Fsts (He > 0.1):

2 populations, normal filtering:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/hist_2pop.jpg" width="400" />

2 populations, random filtering:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/hist_2pop_random.jpg" width="400" />

7 populations, normal filtering:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/hist_7pop.jpg" width="400" />

7 populations, random filtering:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/hist_7pop_random.jpg" width="400" />

Thus, while filtering didn't influence the shape of the Fst distribution, allocating individuals into 2 populations resulted in a very different distribution. This makes sence. In reality, we sampled 7 populations, so that is the real Fst distribution: most sites have Fst near 0.02, low, but not 0 as they are a bit more similar within group then between group. However, when I put the individuals into 2 populations (based on pH) most differences disappear as there are many many sites with Fst = 0.

(Note: there were a total of 991,431 SNPs, for 2 pops 358,803 SNPs were He < 0.1, for 7 pops 357,221 SNPs with He < 0.1, so histograms were made with approximately the same number of SNPs, even if it doesn't look like it)

#### Using OutFlank:

Running OutFlank of default settings (as written above in code chunk) resulted in a relatively good fit for both of the distributions.

OutFlank fit:

2 pops:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/outflank_2pop_default.jpg" width="400" />

7 pops:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/outflank_7pop_default.jpg" width="400" />

Distribution of p-values:

2 pops:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/p_hist_2pop_default.jpg" width="400" />

note the excess of p-values near 1, indicating poor fit of the left tail. Unfortunately, this is a fault of OutFlank, "OutFLANK will not fit the left tail of the FST distribution well".

7 pops:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/p_hist_7pop_dafault.jpg" width="400" />

note that this distribution is uniform, indicating a good fit.

#### Outliers:

For the 7 populations one: no outliers on default settings, and not even with qthreshold=0.1. thus, modified it the following way:

Hmin = 0.05, qthreshold=0.1 ->123 outliers, Hmin = 0.05, qthreshold=0.1, righttrimfraction=0.1 -> 215

Fit: (modified the above two ways leads to same plot)

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/outflank_7pops_q01_h005.jpg" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/p_hist_7pop_q01_h005.jpg" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/outliers_7pop_q01_h005_r01.jpg" width="400" />

For the 2 populations one: 2241 outliers with default settings  (q, Hmin, trims, everything) 

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT/outliers_2pops_default.jpg" width="400" />

I was worried about the distribution of my Fsts, the fit, and the p-value distribution. However, I found this tutorial where they had a similar distribution to mine and they still used it. to be fair the did LD pruning to make the fit better, but I am not gonna do that because LD in urchins is very low (https://rpubs.com/lotterhos/OutFLANK_cichlid_pruning), (LD citation: https://www.biorxiv.org/node/126985.full). Also, another paper used outflank with the same looking distribution with the same kind of fit: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008119.

#### For analysis  (GO enrichment, Chi-squared, etc.) of outliers above and outliers coming from other tests (e.g. LFMM) visit : [this markdown](https://github.com/Cpetak/urchin_adaptation/blob/main/Analysis_of_outliers.md)
