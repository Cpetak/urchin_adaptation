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

### Results - 2 pops

#### Annotating outliers with v5.0

I annotated the 2241 outlier SNPs using the following code: WGS/annotate_outs/annotate.py, process_raw available [here.](https://github.com/Cpetak/urchin_adaptation/blob/main/code/process_raw.py)

```bash
awk -F "," '{print $2}' top_fst_2pops_default.csv > top_fst_2pops_loci # get only relevant column
awk -F "," '$1 != "NA"' top_fst_2pops_loci > fixed_top
mv fixed_top top_fst_2pops_loci
```

```python
import pandas as pd
from collections import Counter

#PROCESSING DATA

#process_raw and data files need to be in the folder you are running this script in

import process_raw #import process_raw.py, also on github

loci = pd.read_csv('top_fst_2pops_loci')
df = pd.read_csv('all_annotations.txt', sep="\t", header=None, names=["chr","source","type","start","stop","s1","s2","s3","info"])
df.drop(df[df['type'] == "region"].index, inplace = True)
df['info'] = df['info'].str.replace(',',';') #will help with parsing region info

# first column of loci df needs to be chr_pos for program to work
#loci["posi"]=loci["posi"].astype(int)
#loci["strpos"]=loci["posi"].astype(str)
#loci["ID"]=loci["chr"]+"_"+loci["strpos"]
#cols = list(loci.columns)
#cols = [cols[-1]] + cols[:-1]
#loci = loci[cols]

#loci["posi"]=loci["posi"].astype(str)

print(loci.head())

#ANNOTATING

annot_df=process_raw.annotate_raw(df,loci) #function defined in process_raw.py, goes through dataframe row by row and find position in annotation file downloaded from ncbi (.gff). available here: https://www.dropbox.com/s/ontctxfee9x7fe4/all_annotations.txt?dl=0
annot_df.pos = annot_df.pos.astype(int)
annot_df.to_csv("annotated01.csv") #basic annotation, includes overlaps
annot_nooverl=annot_df[annot_df['region']!="genes_overlap"]
print("done basic annotation")

overl = annot_df[annot_df['region']=="genes_overlap"]
overlapping_loci=process_raw.process_overlap(df,overl)
cdf=pd.concat([annot_nooverl,overlapping_loci]) #for now we keep both genes if SNPs falls in both
cdf.pos = cdf.pos.astype(int)
cdf.to_csv("annotated02.csv") #basic annotation, overlaps are resolved - note: this results in 2 or more rows in output csv for the same position
print("done processing overlaps annotation")

promoters=process_raw.process_annotation_data(df) #finds positions in promoter regions (within 5000bp of TSS)
missing_annot = annot_df[annot_df['region']=="not_annot"]
promoter_loci=process_raw.process_promoter(missing_annot,promoters)
promoter_loci.pos = promoter_loci.pos.astype(int)
print("done promoter processing")

#CHECK POINT

only_promoter_loci=promoter_loci[promoter_loci['region']!="not_annot"] #loci in promoters
only_not_annot_prom=promoter_loci[promoter_loci['region']=="not_annot"] #loci not annotated even after promoter processing
only_annot_new_df=annot_df[annot_df['region']!="not_annot"] #loci anywhere else

check1=len(only_promoter_loci) + len(only_annot_new_df) + len(only_not_annot_prom)
check2=len(loci)
print("is everything in order?")
print(check1)
print(check2)

only_promoter_loci.to_csv("promoters.csv") #includes all SNP that fall in promoters
print("done done")
```

-> 1294 posi were annotated, 179 fell into promoters

 [annotation for all loci](https://github.com/Cpetak/urchin_adaptation/blob/main/data/annotated02_fst_2pops.csv) (overlaps resolved)

 [promoters](https://github.com/Cpetak/urchin_adaptation/blob/main/data/promoters_fst_2pops.csv)

Converted LOC to SPU using the GenePageGeneralInfo_AllGenes.txt file available on the echinobase.

Code to make a mapping between LOC and SPU based on above file:

```python
import pandas as pd
import json

df=pd.read_csv("SPU_LOC_echino.txt", sep="\t", header=None, names=["ID","LOC","name","gen","test"]) # same as GenePageGeneralInfo_AllGenes.txt

# building dictionary with SPU
SPU_LOC={}
for index,row in df.iterrows():
  newlocs=[]
  newspus=[]
  for coli in ["ID","LOC","name","gen","test"]: #columns where there can be gene name
    if "LOC" in str(row[coli]): #if you find LOC, there could still be SPU
      stuff = row[coli].split("|")
      if len(stuff) > 1:
        locs=[i for i in stuff if "LOC" in i]
        for l in locs:
          newlocs.append(l.split(" ")[-1])
        spus= [i for i in stuff if "SPU" in i]
        for s in spus:
          newspus.append(s)
      else:
        newlocs.append(row[coli].split(" ")[-1])

    if "SPU" in str(row[coli]): #if you find SPU, there could still be LOC
      x = row[coli].split("|")
      locs=[i for i in x if "LOC" in i]
      if len(locs)>0:
        for l in locs:
          newlocs.append(l.split(" ")[-1])
      y=[i for i in x if "SPU" in i]
      for s in y:
        newspus.append(s)

  for l in set(newlocs):
    for s in set(newspus):
      if l in SPU_LOC.keys():
        SPU_LOC[l].append(s)
      else:
        SPU_LOC[l]=[s]

with open('data.json', 'w') as fp:
    json.dump(dict, fp)

```

dictionary is available here: TODO

Then I used this dictionary to map my LOCs to SPUs. Problem: in a lot of cases more than 1 SPU corresponds to a LOC. For Enrichment, I'll use just one SPU at random (choose another if the program doesn't recognise it or doesn't have GO associated to it) - i am trying to see if there is a better way to choose the SPU

TODO code for converting my LOC to SPU, also output file

#### GO enrichment:

1. First, I tried using this tool: http://geneontology.org/, BUT 99% of the SPUs were not recognised, classified as "unclassified"
2. Then, I realised that the legacy echinobase has the GO term for each gene, this made a program to retrieve that. TODO get code here as well as output file -> now I have GO term for each of my SPU
3. Realised that if I don't use an online tool with Spur info already loaded in above info is not enough, all GO term and associated gene information will be needed. oh well, maybe somehow above file will still be useful in the future.
4. Found out that https://string-db.org/ can also do enrichment, all kinds, not just GO. I loaded in all SPUs (including duplicates) -> "there were no significant pathway enrichments observed in the following categories:
   Biological Process (Gene Ontology), Molecular Function (Gene Ontology), Cellular  Component (Gene Ontology), Reference publications (PubMed), KEGG  Pathways, Reactome Pathways, WikiPathways, Protein Domains and Features  (InterPro)." But there were other enrichments, see screenshots on Desktop. Upon further inspection I realised that this is because they have very few genes associated with GO term in their database (a file in the downloads option). Same when I use the LOC gene names...
5. Found this link really useful so check back here if all else fails: https://www.biostars.org/p/261816/
6. Found a python package that looks easy to use, but as I mentioned in 3., I will need all GO terms with all associated genes to use this. COME BACK HERE when I have that, this is my best option right now. Plus comments in the biostars link. https://github.com/atarashansky/LightGOEA
7. I could just wait for them to email me all the info I need in 6. BUT I found the uniprot has that info too! So now I am trying to download that... DOWNLOADED! https://www.ebi.ac.uk/QuickGO/annotations?taxonId=7668&taxonUsage=exact Yey, it seems like it is what i needed. now i just have to go from LOC to uniprotKB... uniprotKB has that information!!! on their gene page there is both uniprot ID and LOC. This is awesome because then I don't have to convert to SPU AND I don't have to use legacy.echinobase AND I don't have to wait for the data to be email to me... I'll still check it tho when it arrives. 

code you might be looking for: https://colab.research.google.com/drive/1EV-LnwYhuEopIt9fyoZSLK9eFK0QJr54?usp=sharing

TODO clean desktop, home and downloads

#### Finding interesting genes:

For searching in lists of biomineralisation, differentially expressed, etc genes, I'll use all SPU alternatives.

- Biomineralisation genes: from Melissa's paper,  https://pubmed.ncbi.nlm.nih.gov/28141889/, IND
- Differentially expressed genes:
  - South vs North after common garden conditions, https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12036 IND
  - Expression of genes of different populations in response to low pH, https://pubmed.ncbi.nlm.nih.gov/28141889/
    - Genes differentially expressed between S. purpuratus populations following one day of exposure to low pH seawater, GENES UP-REGULATED IN POPULATIONS MOST FREQUENTLY EXPOSED TO PH <7.8 & DOWN-REGULATED IN POPULATIONS LESS FREQUENTLY EXPOSED TO PH <7.8, (DE_1.csv), GENES DOWN-REGULATED IN POPULATIONS MOST FREQUENTLY EXPOSED TO  PH <7.8 & UP-REGULATED IN POPULATIONS LESS FREQUENTLY EXPOSED TO PH  <7.8 (DE_2.csv)
    - Genes differentially expressed between S. purpuratus populations following seven days of exposure to low pH seawater, GENES UP-REGULATED IN POPULATIONS MOST FREQUENTLY EXPOSED TO PH  <7.8 & DOWN-REGULATED IN POPULATIONS LESS FREQUENTLY EXPOSED TO PH  <7.8 (DE_3.csv), GENES DOWN-REGULATED IN POPULATIONS MOST FREQUENTLY EXPOSED TO  PH <7.8 & UP-REGULATED IN POPULATIONS LESS FREQUENTLY EXPOSED TO PH  <7.8, (DE_4)
- Genes shown to be related to ph (i.e. SNPs correlated to pH conditions of pops), https://academic.oup.com/icb/article/53/5/857/733987 IND
- Artificial selection low vs normal pH, 1 vs 7 days, allele frequencies that changed, https://www.pnas.org/content/110/17/6937, IND
- Important genes:
- Transcription factors:

HOL TARTOK: I was going through my Melissa_supplementary folder, finished with 2017_expression, next up is 2019_extreme_alleles

ALSO: I need to add these files above to GitHub from the downloads folder where they are currently at

#### Annotating outliers with v3.1

Used this tool to convert locations to v3.1: https://www.ncbi.nlm.nih.gov/genome/tools/remap

Most successfully remapped: 54 failed out of 2241

[remapped loci](https://github.com/Cpetak/urchin_adaptation/blob/main/data/annotated02_fst_2pops_remappedto3.1.txt) 
