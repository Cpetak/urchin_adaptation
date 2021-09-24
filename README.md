# Urchin local adaptation
This repo was made to cleanly demonstrate how I got from raw NGS data to different sets of loci putatively under positive selection using Fst outlier, Bayenv and reduced nucleotide diversity measures.

## Step 1: From reads to bams
<details>
  <summary>Click to view detailed code</summary>

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
covdata <- read.table("covs.txt") #this is without outliers
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
	
## Step 3: Getting per-site Fst values - global, #TODO

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
-minInd 17 
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

