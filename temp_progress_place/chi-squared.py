
"""# Chi-squared"""

from scipy.stats import chisquare
import pandas as pd
import numpy as np

df=pd.read_csv("gathered_annotation.csv")
additional_lnc=pd.read_csv("lncdf.csv")
notannot_enhancers=pd.read_csv("reg_regions_withconfidence.csv")
intronic_enhancers=pd.read_csv("reg_regions_withconfidence_introns.csv")

"""first, the two extremes

## Method 1: overlaps count

Note: redo when I have enhancer info and additional lncRNA
"""

num_exon=len(df[df["exon"]!=0]) # could have overlapped with something else as well...
num_3UTR=len(df[df["3'UTR"]!=0])
num_5UTR=len(df[df["5'UTR"]!=0])
num_intron=len(df[df["intron"]!=0])
num_other=len(df[(df["pseudogene"]!=0) & (df["alternative UTR"]!=0)])
num_lnc=len(df[df["lnc_RNA"]!=0]) + len(additional_lnc)
num_enhancer=len(notannot_enhancers[notannot_enhancers["confidence"]>0]) + len(intronic_enhancers[intronic_enhancers["confidence"]>0])
num_promoter=len(df[df["promoter"]!=0]) # could overlap with other promoter as well as gene
num_noncoding=len(df[df["Not_annot"]==1])
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding] # num_other?

all_all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding, num_other] # num_other?

print((sum(all_all_hit)-len(df))/len(df)*100) # it is this percent of an increase due to my method of counting

#calculated based on https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Strongylocentrotus_purpuratus/102/
#count times mean length (bp)

bp_exon=87323990 # double checked
bp_3UTR=40771628 # double checked
bp_5UTR=10192907 # double checked
bp_intron=412196958 # double checked
bp_lnc_paper= # number of basepairs (no double counting due to overlaps) from paper
bp_lnc_ncbi=4850968 #only from ncbi
bp_promoter=32947 * 2000 #this is based on the hypothesis that promoter = 2kb upstream gene
bp_noncoding=921855793-bp_exon-bp_lnc-bp_promoter-bp_intron-bp_3UTR-bp_5UTR
all_len=[bp_exon,bp_3UTR,bp_5UTR,bp_intron,bp_lnc,bp_promoter,bp_noncoding]

all_percent_expected=[]
for a in all_len:
  all_percent_expected.append(a/921855793) #percentage of genome corresponding to each region

#number of expected loci in each region if the high Fst loci were randomly distributed between the regions
all_expected=[round(element * sum(all_hit)) for element in all_percent_expected]

print(chisquare(all_hit, f_exp=all_expected))
print("[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding]")
print(all_hit)
print(all_expected)
print([i / j for i, j in zip(all_hit, all_expected)])

"""## Method 2: everything is counted once, more conservative

gene-gene neither is counted, promoter-promoter counted as 1, gene-promoter and gene-promoter-promoter only gene is counted.
"""

# Make alternative df
df2 = pd.read_csv("annotated01.csv") #with overlaps in 1 row, called genes_overlap
df2 = df2.iloc[: , 1:]
df2 = df2[["chr", "pos", "region","gene"]]
print(len(df2[df2["region"]=="genes_overlap"]))
df2 = df2[df2["region"]!="genes_overlap"]

p2=pd.read_csv("promoters02.csv")
p2 = p2.iloc[: , 1:]

# Not annotated
not_annot=p2[p2["region"]=="not_annot"]
df2=df2.merge(not_annot, how='outer', on=['chr', 'pos'])
df2.drop('gene_y', axis=1, inplace=True)
df2=df2.rename(columns={"region_y": "Not_annot"})
df2.Not_annot = df2.Not_annot.fillna('NaN')
df2["Not_annot"] = df2["Not_annot"].map({"not_annot":1, "NaN":0,})

# Promoters
proms=p2[p2["region"]!="not_annot"]
proms=proms[proms["region"]!="genes_overlap"]
proms["region"] = proms["region"].map({"protein_coding":"protein_coding", "lnc_RNA":"lnc_RNA","tRNA":"tRNA","intron":"protein_coding" })
proms2=proms.groupby(["chr","pos"]).agg(lambda x: list(x))

df2=df2.merge(proms2, how='outer', on=['chr', 'pos'])
df2=df2.rename(columns={"gene": "prom_gene", "region":"prom_region"})
#df['promoter'] = df["promoter"].where(df["promoter"].isnull(), 1).fillna(0).astype(int)

df2["promoter"]=df2["prom_region"]
df2['promoter'] = df2['promoter'].str.len()
df2['promoter'] = df2['promoter'].fillna(0)

df2

num_exon=len(df2[df2["region_x"]=="exon"])
num_3UTR=len(df2[df2["region_x"]=="3'UTR"])
num_5UTR=len(df2[df2["region_x"]=="5'UTR"])
num_intron=len(df2[df2["region_x"]=="intron"])
num_other=len(df2[(df2["region_x"]=="pseudogene") & (df2["region_x"]=="alternative UTR")])
num_lnc=len(df2[df2["region_x"]=="lnc_RNA"]) #TODO add extra lncRNAs here
# missing rRNA, tRNA etc
num_promoter=len(df2[df2["promoter"]!=0]) # could overlap with other promoter as well as gene
num_noncoding=len(df2[df2["Not_annot"]==1])
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding] # num_other?

all_all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding, num_other] # num_other?

print((sum(all_all_hit)-len(df2))/len(df2)*100) # it is this percent of an increase due to my method of counting

#calculated based on https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Strongylocentrotus_purpuratus/102/
#count times mean length (bp)

bp_exon=87323990
bp_3UTR=40771628
bp_5UTR=10192907
bp_intron=412196958
#num_other??
#bp_lnc=15513236 #(including both annotation data from ncbi AND extra lnc_RNA file)
bp_lnc=4850968 #only from ncbi
bp_promoter=32947 * 2000
bp_noncoding=921855793-bp_exon-bp_lnc-bp_promoter-bp_intron-bp_3UTR-bp_5UTR
all_len=[bp_exon,bp_3UTR,bp_5UTR,bp_intron,bp_lnc,bp_promoter,bp_noncoding]

all_percent_expected=[]
for a in all_len:
  all_percent_expected.append(a/921855793) #percentage of genome corresponding to each region

#number of expected loci in each region if the high Fst loci were randomly distributed between the regions
all_expected=[round(element * sum(all_hit)) for element in all_percent_expected]

print(chisquare(all_hit, f_exp=all_expected))
print("[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding]")
print(all_hit)
print(all_expected)
print([i / j for i, j in zip(all_hit, all_expected)])

"""## Method 3: same as method 2 but for gene-gene overlap select 1 at random instead of discarding all

Technically it drops first duplicate but since they were added to the csv randomly it shouldn't matter, so practically we end up selecting randomly
"""

# Make alternative df
df2 = pd.read_csv("annotated02.csv")
df2 = df2.iloc[: , 1:]
df2 = df2[["chr", "pos", "region","gene"]]
df2=df2.drop_duplicates(subset=['chr','pos']) # keep only one for overlaps

p2=pd.read_csv("promoters02.csv")
p2 = p2.iloc[: , 1:]

# Not annotated
not_annot=p2[p2["region"]=="not_annot"]
df2=df2.merge(not_annot, how='outer', on=['chr', 'pos'])
df2.drop('gene_y', axis=1, inplace=True)
df2=df2.rename(columns={"region_y": "Not_annot"})
df2.Not_annot = df2.Not_annot.fillna('NaN')
df2["Not_annot"] = df2["Not_annot"].map({"not_annot":1, "NaN":0,})

# Promoters
proms=p2[p2["region"]!="not_annot"]
proms=proms[proms["region"]!="genes_overlap"]
proms["region"] = proms["region"].map({"protein_coding":"protein_coding", "lnc_RNA":"lnc_RNA","tRNA":"tRNA","intron":"protein_coding" })
proms2=proms.groupby(["chr","pos"]).agg(lambda x: list(x))

df2=df2.merge(proms2, how='outer', on=['chr', 'pos'])
df2=df2.rename(columns={"gene": "prom_gene", "region":"prom_region"})
#df['promoter'] = df["promoter"].where(df["promoter"].isnull(), 1).fillna(0).astype(int)

df2["promoter"]=df2["prom_region"]
df2['promoter'] = df2['promoter'].str.len()
df2['promoter'] = df2['promoter'].fillna(0)

num_exon=len(df2[df2["region_x"]=="exon"])
num_3UTR=len(df2[df2["region_x"]=="3'UTR"])
num_5UTR=len(df2[df2["region_x"]=="5'UTR"])
num_intron=len(df2[df2["region_x"]=="intron"])
num_other=len(df2[(df2["region_x"]=="pseudogene") & (df2["region_x"]=="alternative UTR")])
num_lnc=len(df2[df2["region_x"]=="lnc_RNA"]) #TODO add extra lncRNAs here
# missing rRNA, tRNA etc
num_promoter=len(df2[df2["promoter"]!=0]) # could overlap with other promoter as well as gene
num_noncoding=len(df2[df2["Not_annot"]==1])
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding] # num_other?

all_all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding, num_other] # num_other?

print((sum(all_all_hit)-len(df2))/len(df2)*100) # it is this percent of an increase due to my method of counting

#calculated based on https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Strongylocentrotus_purpuratus/102/
#count times mean length (bp)

bp_exon=87323990
bp_3UTR=40771628
bp_5UTR=10192907
bp_intron=412196958
#num_other??
#bp_lnc=15513236 #(including both annotation data from ncbi AND extra lnc_RNA file)
bp_lnc=4850968 #only from ncbi
bp_promoter=32947 * 2000
bp_noncoding=921855793-bp_exon-bp_lnc-bp_promoter-bp_intron-bp_3UTR-bp_5UTR
all_len=[bp_exon,bp_3UTR,bp_5UTR,bp_intron,bp_lnc,bp_promoter,bp_noncoding]

all_percent_expected=[]
for a in all_len:
  all_percent_expected.append(a/921855793) #percentage of genome corresponding to each region

#number of expected loci in each region if the high Fst loci were randomly distributed between the regions
all_expected=[round(element * sum(all_hit)) for element in all_percent_expected]

print(chisquare(all_hit, f_exp=all_expected))
print("[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_promoter,num_noncoding]")
print(all_hit)
print(all_expected)
print([i / j for i, j in zip(all_hit, all_expected)])
