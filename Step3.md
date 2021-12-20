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

## Use Outflank for per-site Fst

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
awk -F, '$3 > 0.1' fixed_combined.csv > fixed_combined_goodhe.csv #getting only sites with He > 0.1
cat fixed_combined_goodhe.csv | cut -d ',' -f2,7 > twocol.csv #keeping only position and FSTNoCorr columns
sort -k 2 -t , -n -r twocol.csv > sorted_twocol.csv #sort by Fst
cat sorted_twocol.csv | grep -v "e" > test.csv remove first few lines that were incorrectly sorted due to eg e-10
head -6823 test.csv > reduced_sorted_twocol.csv #keep only high Fst pos, the threshold will depend on the distribution
```

Where Pop.txt 

