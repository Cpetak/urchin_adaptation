# Urchin local adaptation
This repo was made to cleanly demonstrate how I got from raw NGS data to different sets of loci putatively under positive selection using Fst outlier, Bayenv and reduced nucleotide diversity measures.

## Step 1: From reads to bams

[Code for Step 1](https://github.com/Cpetak/urchin_adaptation/blob/main/Step1.md)

## Step 2: from bams to genotype likelihoods and PCA
[Code for Step 2](https://github.com/Cpetak/urchin_adaptation/blob/main/Step2.md)

Getting the filtered sites

[Code for Step 2](https://github.com/Cpetak/urchin_adaptation/blob/main/Step2.md)



## Step 3: Getting per-site Fst values - global
[Code for Step 3](https://github.com/Cpetak/urchin_adaptation/blob/main/Step3.md)


## Step 4: Running LFMM

[Code for Step 4](https://github.com/Cpetak/urchin_adaptation/blob/main/Step4.md)

## Step 5: PCAngsd

[Code for Step 5](https://github.com/Cpetak/urchin_adaptation/blob/main/Step5.md)

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

