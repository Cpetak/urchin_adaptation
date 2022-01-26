# Analysis of outliers

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



8. I will use topGo and I can convert LOC to UniprotID using this tool: https://www.uniprot.org/uploadlists/. select From Ensemble Genomes To uniprot

9. used this code to transform uniprot - GO output file into the mapping file topGO expects:

   ```bash
   awk -F "\t" '{print $2"\t"$5}' QuickGO-annotations-1642716310981-20220120.tsv > temp_mapping
   sed '$!N; /^\(.*\)\n\1$/!P; D' temp_mapping > temp2_mapping # It deletes duplicate, consecutive lines from a file
   awk 'BEGIN{FS="\t"} {for(i=2; i<=NF; i++) { if (!a[$1]) a[$1]=$1FS$i ;else a[$1]=a[$1]","$i};if ($1 != old) b[j++] = a[old];old=$1 } END{for (i=0; i<j; i++) print b[i] }' temp2_mapping > GO_mapping_topGO #it collapses repeated lines into 1, comma separated
   ```

10. translated LOCs to uniprot IDs using above link

"934 out of 956 Ensembl Genomes identifiers were successfully mapped to 1769 UniProtKB  IDs in the table below." See colab for how i made the uniprotID list of "intereting genes"



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