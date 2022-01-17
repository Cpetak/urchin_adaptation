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

For all figures, see [this folder](https://htmlpreview.github.io/?https://github.com/Cpetak/urchin_adaptation/blob/main/images/FCT) 



however, like this, fit is poor and p value distribution is not uniform

so trying with left and right trim fraction parameters -> nothing changed

So I created a random subset of vcf_tail and rerun everything above on that to see if my filtering is biasing the distribution (not chi-squared) 

-> same distribution of Fst, now trying with 7 pops instead of 2 -> different distribution? outflank works with this?

I might end up just taking the top 1% as before... with bootstrapping

amúgy elbasztam előbb úgyhogy újra kell futtatni az outflank-et, k=2 nem k=7
