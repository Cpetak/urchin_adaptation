# Analysis of outliers

## 2 pops per-site Fst values from Step 3

### Chi-squared for regions

#### Annotating outliers with v5.0 for nonregulatory variation

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

#ANNOTATING

annot_df=process_raw.annotate_raw(df,loci, one_length_print=False) #function defined in process_raw.py
annot_df.pos = annot_df.pos.astype(int)
annot_df.to_csv("annotated01.csv") #basic annotation, includes overlaps in one row
annot_nooverl=annot_df[annot_df['region']!="genes_overlap"]
print("done basic annotation")

overl = annot_df[annot_df['region']=="genes_overlap"]
overlapping_loci=process_raw.process_overlap(df,overl)
cdf=pd.concat([annot_nooverl,overlapping_loci]) #for now we keep both genes if SNPs falls in both
cdf.pos = cdf.pos.astype(int)
cdf.to_csv("annotated02.csv") #basic annotation, overlaps are resolved into multiple rows
print("done processing overlaps annotation")

promoters=process_raw.process_annotation_data(df)
missing_annot = annot_df[annot_df['region']=="not_annot"]
promoter_loci=process_raw.process_promoter(missing_annot,promoters)
promoter_loci.pos = promoter_loci.pos.astype(int) # loci in promoters, includes promoter-promoter overlaps in one row
print("done promoter processing")

poverl = promoter_loci[promoter_loci['region']=="genes_overlap"]
poverlapping_loci=process_raw.process_overlap(promoters,poverl)
cdf2=pd.concat([promoter_loci,poverlapping_loci])
cdf2.pos = cdf2.pos.astype(int)
cdf2.to_csv("promoters02.csv") # loci in promoters, promoter-promoter overlaps are resolved into mutliple rows
print("done processing overlaps in promoter annotation")

# we accounted for pos that fall in more than 1 gene region, and pos that fall in more than 1 promoter region
# now we correct for if pos falls into both gene body and promoter region

new_df=process_raw.annotate_raw(promoters,loci, one_length_print=True) # redo annotation with with +5000 for promoters
new_df = new_df[new_df['region']=="genes_overlap"] #now it can be gene-gene, gene-promoter, gene-gene-promoter, gene-promoter-promoter or promoter-promoter overlap
new_df=new_df[~new_df.pos.isin(poverl.pos)] #remove promoter-promoter overlap
new_loci=process_raw.process_overlap(promoters,new_df)

ori_overlaps=annot_df[annot_df["region"]=="genes_overlap"]
oridf=cdf[cdf.pos.isin(ori_overlaps.pos)] # get annotations for overlapped genes

new_loci.gene = new_loci.gene.fillna('NaN')
oridf.gene = oridf.gene.fillna('NaN')
annot_df.gene = annot_df.gene.fillna('NaN')

new_loci.to_csv("new_loci.csv")
print(new_loci.head())

print("starting gene-promoter overlap resolve")
outdf=pd.DataFrame()
for index, row in new_loci.iterrows():
  found_in_oridf = oridf[(oridf["pos"]==row[2]) & (oridf["chr"]==row[0])] # chr, pos also in origianl df containing overlaps (no promoter)
  others_in_new_loci = new_loci[(new_loci["pos"]==row[2]) & (new_loci["chr"]==row[0])] # other same chr, pos in this df
  if len(found_in_oridf) > 0: # gene-gene overlap?
    if len(found_in_oridf) == len(others_in_new_loci): # regular gene-gene overlap, already in annotated02.csv
      pass
    elif len(found_in_oridf) < len(others_in_new_loci): #special case, gene-gene-promoter, resolve later if it comes up
      print("Error: more hits in new_loci than in oridf, due to gene-gene-promoter overlap")
    else:
      print("Error: more hits in oridf than in new_loci")
  else: # gene-promoter or gene-promoter-promoter, this is what we are really interested in
    dfmatch=annot_df[(annot_df["chr"]==row[0]) & (annot_df["pos"]==row[2]) & (annot_df["gene"]==row[1])] # chr, pos, gene name in original df
    if len(dfmatch) == 1: # this is the original gene name, hit for the gene body
      outdf = outdf.append({'chr' : row[0], 'pos' : row[2], 'region' : list(dfmatch.region)[0], "gene" : list(dfmatch.gene)[0]}, ignore_index = True)
    if len(dfmatch) == 0: # this is the new gene hit, hit for the promoter
        r="promoter" + ":" + row[3] # region it would fall in in "extended gene body", if intron, UTR, etc. if is promoter for regular gene
        outdf = outdf.append({'chr' : row[0], 'pos' : row[2], 'region' : r, "gene" : row[1]}, ignore_index = True)
    if len(dfmatch) > 1:
      print("Error: more than one hit in original df")

outdf[['promoter','p_type']] = outdf['region'].str.split(':',expand=True)

print("checking gene-promoter output df")
tdf=outdf.groupby(['chr','pos']).count()
tdf["ok"]= tdf["region"]-tdf["p_type"] # if at least 1 -> there is a nonpromoter region
print(len(tdf[tdf["ok"]<1])) # should be 0, if not, nonpromoter is missing somewhere
print(len(tdf[tdf["p_type"]<1])) # should be 0, if not, promoter is missing somewhere

p_overl_g=outdf[outdf["promoter"]=="promoter"]
p_overl_g["p_type"] = p_overl_g["p_type"].map({"intron":"protein_coding", "3'UTR":"protein_coding", 
                                               "5'UTR":"protein_coding", "alternative UTR":"protein_coding", 
                                               "exon":"protein_coding", "tRNA":"tRNA", "lnc_RNA":"lnc_RNA", "miRNA":"miRNA", "pseudogene":"pseudogene",
                                               "snRNA":"snRNA", "snoRNA":"snoRNA", "rRNA":"rRNA"})
p_overl_g.drop('region', axis=1, inplace=True)
p_overl_g.drop('promoter', axis=1, inplace=True)
p_overl_g["region"]=p_overl_g["p_type"]
p_overl_g.drop('p_type', axis=1, inplace=True)
p_overl_g.to_csv("promoter_gene_overl.csv") # new promoter regions I didn't find before because they are gene-promoter-promoter or gene-promoter overlaps

#CHECK POINT
print("starting check points")
only_promoter_loci=promoter_loci[promoter_loci['region']!="not_annot"] #loci in promoters
only_not_annot_prom=promoter_loci[promoter_loci['region']=="not_annot"] #loci not annotated even after promoter processing
only_annot_new_df=annot_df[annot_df['region']!="not_annot"] #loci anywhere else

check1=len(only_promoter_loci) + len(only_annot_new_df) + len(only_not_annot_prom)
check2=len(loci)
print("is everything in order?")
print(check1)
print(check2)

only_promoter_loci.to_csv("promoters01.csv") #includes all SNP that fall in promoters
print("done done")
```

There are 6 types of annotations:

1. Gene body hits -> annotated01.csv, 
2. Hits in gene body - gene body overlaps are resolved into separate lines -> annotated02.csv
3. Promoter hits -> promoters01.csv
4. Hits in promoter - promoter overlaps resolved into separate lines ->  promoters02.csv
5. Hits in gene body - promoter overlaps  and hits in gene body - promoter - promoter overlaps -> promoter_gene_overl.csv
7. Hits in gene body - gene body - promoter overlaps (rare, this annotation still shows pos as simple gene body - gene body overlap)

 [annotation for all loci](https://github.com/Cpetak/urchin_adaptation/blob/main/data/2_pop_fst_step3/annotated02_fst_2pops.csv) (overlaps resolved) TODO update these links

 [promoters](https://github.com/Cpetak/urchin_adaptation/blob/main/data/2_pop_fst_step3/promoters_fst_2pops.csv)

TODO rerun as I changed promoter = 5000kb to 2000kb

I gathered all of the above annotation data into one big dataframe using the gather_annot_info.py code. This created a dataframe with all outlier loci as rows and columns containing region information such as intron, promoter, etc with a number representing how many of these regions each loci "hit". So for example if a locus was found in 2 overlapping promoter regions, the whole row has 0s except for the promoter column which says 2. From this dataframe, I then extracted rows that were not annotated at all. -> input_for_extra_regannot_notannot.csv, Loci in this file were converted to v3.1 as follows. 

#### Annotating outliers with v3.1 for regulatory regions

Used this tool to convert locations to v3.1: https://www.ncbi.nlm.nih.gov/genome/tools/remap

Opened input_for_extra_regannot_notannot.csv in Visual Studio Code, and from there I copied into space in link above. Had to remove qutation marks and ".1" from the end of the chromosome numbers. Chr, pos space separated.

Most successfully remapped: 54 failed out of 2241, some mapped to multiple regions. #TODO update this

[remapped loci](https://github.com/Cpetak/urchin_adaptation/blob/main/data/annotated02_fst_2pops_remappedto3.1.txt) #TODO update all links in this markdown

Regulatory regions: #TODO find citations for each of these

* ATAC-seq DONE
* DNA-seq DONE
* Chip-seq DONE
* Lvar similarity DONE
* Also Arenas-Mena 2021 DONE https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07936-0
* Enhancer RNA dataset Khor 2021 DONE

Additional: lncRNA DONE

TODO when available: echinobase CRE experimentally validate data, supp of Khor 2021, Computational ID, TFBS

For the above, some files contained overlaps within regions, I resolved these using the check_overlap_within_df.py.

Then, the files above along with the input_for_extra_regannot_notannot_3.1.csv file were the input to the annot_reg_regions.py program. -> output list of positions that fell in an lncRNA region, and another list of positions that fell in any of the regulatory regions instead (ATAC, Chip, L.var, etc.).

#### Putting it together for final analysis

Chi-squared analysis was run as shown in chi-squared.py. It takes gathered_annotation.csv (output of gather_annot_info.py), lncdf.csv and reg_regions_withconfidence.csv (output of annot_reg_regions.py). 

To get the total number of lncRNA nucleotides (for expected number of hits in chi-squared): overlap processed with check_overlap_within_df.py, then gff annotation (processed with filtering_annotation_forlist.py to account for promoters) was subtracted with take_out_gff_from_lncrna.py. Note: sp4.lncRNAs_overlap_processed.csv was first converted from 3.1 to 5.0. Then, each region's length was calculated and summed, together with the number of nucleotides for lncRNA in ncbi. Final number: 17,703,544

To get the total number of enhancer nucleotides (for expected number of hits in chi-squared): overlapped regions from all different resources (ATAC, Chip, etc) into 1 file using the get_enhancer_overlaps.py code (regulatory_regions_lit_review folder, output overlapped_enhancers.csv). Then, I subtracted regions in lncRNA (additional file) and gff (processed as above). Gff file was translated into version 3.1. USED CODE X.

### GO enrichment

1. First, I tried using this tool: http://geneontology.org/, BUT 99% of the SPUs were not recognised, classified as "unclassified"
2. Then, I realised that the **legacy echinobase** has the GO term for each gene, this made a program to retrieve that. https://colab.research.google.com/drive/1dD9OVpyoZBj-E1cDitCeo0_kvJMX3c2e?usp=sharing. Didn't end up running as a) I found a better solution, b) I realsied I need the GO terms for each gene in the urchin genome, not just my genes of interest, for the GO enrichment to work.
3. Found out that https://string-db.org/ can also do enrichment, all kinds, not just GO. TODO. Upon further inspection I realised that this might not work because the string database has very few genes associated with GO term (based on a file in the downloads option). 
4. Found this link really useful so check back here if all else fails: https://www.biostars.org/p/261816/, also this is a potential package to use: https://github.com/atarashansky/LightGOEA
5. To get the GO terms info I need, I emailed people at echinobase but they sent me a file that only had around 300 gene with associated GO terms.
6. It turns out that **Uniprot has GO terms** associated to each gene in the urchin genome, here is the link to retrieve that info: https://www.ebi.ac.uk/QuickGO/annotations?taxonId=7668&taxonUsage=exact. Click on the export button. Saved as QuickGO-annotations-1642716310981-20220120.tsv in this git repo.
7. Used this code to transform uniprot - GO output file into the mapping file topGO expects:

```bash
awk -F "\t" '{print $2"\t"$5}' QuickGO-annotations-1642716310981-20220120.tsv > temp_mapping
sed '$!N; /^\(.*\)\n\1$/!P; D' temp_mapping > temp2_mapping # It deletes duplicate, consecutive lines from a file
awk 'BEGIN{FS="\t"} {for(i=2; i<=NF; i++) { if (!a[$1]) a[$1]=$1FS$i ;else a[$1]=a[$1]","$i};if ($1 != old) b[j++] = a[old];old=$1 } END{for (i=0; i<j; i++) print b[i] }' temp2_mapping > GO_mapping_topGO #it collapses repeated lines into 1, comma separated, output file is in this git repo
```

8. Now that I have all genes with associated GO terms, I need to map **LOC IDs into UniprotIDs** to get my list of interesting genes I can use in the GO enrichment together with the gene to GO mapping file above:
   1. To get the list of LOC: TODO
   2. I can convert LOC to UniprotID using this tool: https://www.uniprot.org/uploadlists/. select From Ensemble Genomes To uniprot.
   3. "934 out of 956 Ensembl Genomes identifiers were successfully mapped to 1769 UniProtKB  IDs in the table below." 
   4. Download tab separated file using the download button -> LOC2uniprot_long.tab TODO
   5. Use code in this colab: https://colab.research.google.com/drive/1lLUHlDrIYj6RE8-ozVoR9PioYPa33FWB?usp=sharing to resolve cases where one LOC mapped to multiple UniprotIDs. Since most UniprotIDs sharing a single LOC mapped to the same GO terms, I just ended up selecting the the UniprotID with the most GO terms. -> list of uniprotIDs of interest
9. Using the the two key files above, I run topgo using the following code:

```R
library(topGO)

geneID2GO <- readMappings("GO_mapping_topGO") # uniprot to GO mapping
geneNames <- names(geneID2GO)

myInterestingGenes <- read.csv("uniprotIDs_fst_2pops.txt", header = FALSE) # list of interesting genes, output of LOC to uniprot mapping
intgenes <- myInterestingGenes[, "V1"]
geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
names(geneList) <- geneNames # geneList but annotated with the gene names

GOdata <- new("topGOdata", 
              ontology = "BP", # ontology of interest (BP, MF or CC)
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") # these are the options I'll be using! checked!
resultFisher

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10) # top 10 enriched terms
allRes

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all') 

```

Echinobase was down but they emailed them and they sent me a file with around 300 genes with associated GO terms, called GeneGoTerms.txt. The genes were in a GENEPAGE number format, so first I translated my LOC list into a GENEPAGE list using the SPU_LOC_echino.txt file and this code: https://colab.research.google.com/drive/19pi5vPVcZvd_OfEzFTj0XaQnnCVox0yc?usp=sharing -> made a json file to map from LOC to GENEPAGE, LOC2GENEPAGE.json, and mapped all my LOCs into GENEPAGE -> GENEPAGE_fst_2pops.txt. Also made a file with all of the GO terms and associated GENEPAGE numbers in a format topGO expects (i.e. the gene universe) -> GO_mapping_topGO_echinoGENEPAGE. This, however, didn't go anywhere as I said, this only included a total of ~300 genes... all of these files are in Echinobase_GO_stuff

### SPU for supplementary data analysis

Converted LOC (the ones with outliers in their gene bodies, no lncRNA but also promoter regions) to SPU using the GenePageGeneralInfo_AllGenes.txt file available on the echinobase.

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

