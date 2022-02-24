# Step 6: Getting per-site Fst values - pair-wise

Folders:

WGS/angsd_new/ for angsd outputs, WGS/angsd_new/pairwise_fst_cleaned for all other files and scripts

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

For FOG and CAP, I also run the above code without outliers -> angsd_noout

Since I have 7 populations, I have 21 possible pairs of populations. For each of the possible pairs:

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

For all pairs containing CAP and FOG the angsd output with no outliers was used (from angsd_noout)

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

Distribution of final pairwise Fst values:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/pairwise_FST_histogram.jpg" width="400" />

Then bootstrapped the top 1% cutoff point using the following code:

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv("subbed.csv")

fsts=df.final_fst.to_numpy()

def bootstrap_sample(amounts):
    return np.random.choice(amounts, len(amounts), replace=True)

def percentile_99(sample):
     return np.percentile(sample, 99)

def bootstrap_confidence_interval(data):
    """
    Creates list of 10000 99th percentile bootstrap replicates.
    """

    itnum=10000

    bs_samples = np.empty(itnum)

    for i in range(itnum):
        bs_samples[i] = percentile_99(bootstrap_sample(data))

    return bs_samples

transactions_ci = bootstrap_confidence_interval(fsts)

print(np.percentile(transactions_ci, 95))

np.savetxt('99s.out', transactions_ci, delimiter=',')

fig3 = plt.figure()
plt.hist(transactions_ci, bins=40)
fig3.savefig('boots.png')
```

Output: 0.0260773276575

```bash
awk -F "," ' $6 >= 0.0260773276575 ' subbed.csv > pair_fst_outs
```

Number of outliers here: 9780

For analysis  (GO enrichment, Chi-squared, etc.) of outliers above and outliers coming from other tests (e.g. LFMM) visit : [this markdown](https://github.com/Cpetak/urchin_adaptation/blob/main/Analysis_of_outliers.md)

