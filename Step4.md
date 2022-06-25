---
typora-copy-images-to: ./images
---

# Step 4: Running LFMM

folders:

WGS/make_vcf/using_vcf/LFMM

converted filtered_vcf (as described in [Code for Step 3](https://github.com/Cpetak/urchin_adaptation/blob/main/Step3.md), section "Filtering vcf_tail.vcf") into LFMM format after removing CAP and FOG outliers:

```bash
# take out first 10 columns - info about position, and columns 31, 57, 58 which are the outlier individuals
cut --complement -d$'\t' -f1,2,3,4,5,6,7,8,9,30,56,57 filtered_vcf > cut_filtered_vcf
# since file is very big I transposed it using python -> num rows=137, num cols=991,430
with open("cut_filtered_vcf") as f:
    raw = f.readlines()
    data = [el.strip().replace("\n","").split("\t") for el in raw]

def get_col(i,data):
    return [el[i] for el in data]

with open("cut_filtered_vcf_t","w") as f:
    for i in range(len(data[0])):
        f.write(" ".join(get_col(i,data))+"\n")
# make into LFMM format
sed 's/0\/0/0/g' cut_filtered_vcf_t | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./9/g' > lfmm_input # 's/NA/9/g' for non-pruned
```

Running LFMM:

Getting K, spack load r@3.6.3

```R
#devtools::install_github("bcm-uga/lfmm")
library(lfmm)
library(tseries)

mydata<-read.matrix("lfmm_input", header = FALSE, sep = " ", skip = 0)

#test<-c(0,0,0,2,0,2,2)
#test_pos<-rep(test, each=20)
#mydata2 <- cbind(mydata, test_pos)

# get genotype data, columns are positions, rows are individuals
#Y <- example.data$genotype
Y <- mydata

#principal component analysis (PCA) can reveal some ‘structure’ in the genotypic data. 
#We perfom PCA by using the prcomp function as follows.
pc <- prcomp(Y)
pdf(file = "PCA.pdf")
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
dev.off()
#For the A. thaliana samples, the screeplot indicates that there are around K=6 main components in the data. 
#We will use K=6 latent factors in subsequent analyses of the A. thaliana example.
```



Getting p-values

```R
library(lfmm)
library(tseries)
Y<-read.matrix("lfmm_input", header = FALSE, sep = " ", skip = 0)
#order of vcf: BOD,CAP,FOG,KIB,LOM,SAN,TER
BOD <- rep(0.05026657, each=20)
CAP <- rep(0.07188418, each=19)
FOG <- rep(0.1780976, each=18)
rest <- rep(c(0.003375068,0.00729927,0,0), each=20)
X <- c(BOD, CAP, FOG, rest)

#The R package contains two main functions for estimating the parameters of LFMMs: 
# ridge_lfmm and lasso_lfmm. 
# ridge minimizes sum of squared residuals + lambda * slope squared -> less sensitive to changes in x, lambda is a scaler, the larger lambda the less sensitive y is to x (i.e. flatter the fitted line)
# overall shrinks parameters, making our predictions less sensitive to them
# vs lasso regression: sum of squared residuals + lambda * absolute_value(slope), again making prediction (y) less sensitive to x
# lasso can shrink parameters to 0, ridge can't (when increasing lambda) -> lasso can exclude useless variables from our model
# ridge is more useful when most parameters are useful -> are most of my parameters useful? I think so!
# so I need ridge... let's do both tho and compare results?
# ridge is also the one used in most papers citing the lfmm paper
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 7) #determined from PCA, TODO

#The ridge_lfmm function returns an object that contains the latent variable score matrix U, 
# the latent variable loading matrix U, and B the effect sizes for all SNPs.

## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif") #gif is default, and what others use
# above funtion's documentation: https://rdrr.io/github/bcm-uga/lfmm/man/lfmm_test.html. "Additional corrections are required for multiple testing"

pvalues <- pv$calibrated.pvalue 
write.csv(pvalues,"LFMM_ridge_pvalues.csv")
print("created csv with p-vals")

#The histogram of test significance values is expected to be flat, with a peak near zero. 
#A QQ-plot is displayed as follows. Explanation of QQ-plot: https://towardsdatascience.com/q-q-plots-explained-5aa8495426c0

png(file = "QQplot.png")
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
dev.off()

## Manhattan plot
png(file = "Manhattan.png")
plot(-log10(pvalues), 
     pch = 19, 
     cex = .9, 
     xlab = "SNP", ylab = "-Log P",
     col = "black")
dev.off()
```

Tried above code with k=1,2,7, and both binary and continous variables for pH.

Then, run the following code for correction of multiple testing:

```R
library(qvalue)

results<-read.csv(file="k_7_cont/LFMM_ridge_pvalues.csv", header=TRUE)
print(head(results))

pvalues<-results$V1

qobj <- qvalue(p = pvalues)
qvalues <- qobj$qvalues
pi0 <- qobj$pi0
lfdr <- qobj$lfdr

summary(qobj)

pdf(file = "qvals.pdf")
hist(qobj)
dev.off()
```

k=7, binary -> no significant qval

k=7, continous -> 30 under 0.05

k=1, binary ->  no significant qval

k=2, binary -> no significant qval MAKES SENSE THAT CONTINOUS IS BETTER!

k=2, continous -> 7 under 0.05, 29 under 0.1, for pval 7681 under 0.01, 1242 under 0.001

k=1, continous -> for qval 8 under 0.05, 25 under 0.1, for pval 7490 under 0.01, 1149 under 0.001

so I am choosing continous, and k=1 or k=2 doesn't matter, p-values are very similary, qqplots nearly identical. No k=7 as it doesn't make sence based on PCA plot and also qqplot looks off

## Fix idea number 1:

In WGS/make_vcf/using_vcf/LFMM_pruned

LD pruning. [Code for LD pruning here](https://github.com/Cpetak/urchin_adaptation/blob/main/LD pruning.md) 

0.5 -> 215,911 of 991,430 variants removed - 20%

0.2 -> 393,874 of 991,430 variants removed - 40%

neither changes LFMM results, nor did normalising env data

Tested LFMM by adding a position with prefect allele frequencies and it gave me a p value of 0 basically and found it as significant by qvalues so that is not the issue… test01: StoN: 0,0,0,50%,0,50%,100%, test02: StoN: 0,0,0,100%,0,100%,100%

## Fix idea number 2:

In WGS/make_vcf/using_vcf/LFMM_testing_env

I am not using good values to represent pH variability? I reanalised the pH data and came up with other measures to correlate allele frequencies to. Code is [here](https://colab.research.google.com/drive/1rIY1rALjo2Rom6GuMkxzH9pBpk8xXEWT?usp=sharing), it is a colab notebook. These are:

No date restrict means all data was considered in the calculation, date restrict means that only data recorded in April, May and June are considered, as well as only years 2011 and 2012. This second option is reasonable to the following reasons: 1) that is the period for larval life stage, 2) least missing data, 3) most even distribution of datapoints between populations.

|      | no date restrict  | BOD      | CAP      | FOG      | KIB      | LOM      | SAN (spec) | TER     |
| ---- | ----------------- | -------- | -------- | -------- | -------- | -------- | ---------- | ------- |
| 1    | freq7.8           | 0.050751 | 0.071626 | 0.166984 | 0.003593 | 0.007647 | 0          | 0       |
| 2    | mins              | 7.539    | 7.555    | 7.541    | 7.73     | 7.767    | 8          | 7.819   |
| 3    | ave100            | 7.59296  | 7.62761  | 7.56631  | 7.8011   | 7.82439  | 8          | 7.84115 |
| 4    | lower1%           | 7.676    | 7.677    | 7.631    | 7.82384  | 7.80476  | 8          | 7.864   |
|      | **date restrict** |          |          |          |          |          |            |         |
| 5    | freq7.8           | 0.001732 | 0.040543 | 0.008983 | 0.005094 | 0.002165 | 0          | 0       |
| 6    | mins              | 7.775    | 7.597    | 7.668    | 7.73     | 7.795    | 8          | 7.819   |
| 7    | ave100            | 7.81032  | 7.70106  | 7.7493   | 7.82525  | 7.84692  | 8          | 7.84115 |
| 8    | lower1%           | 7.82791  | 7.71876  | 7.807    | 7.829    | 7.80885  | 8          | 7.864   |

- Freq7.8 = number of datapoints under 7.8, standardized by the total number of datapoints recorded at that location.
- Mins = minimum pH value recorded at that location.

- Ave100 = lowest 100 pH values averaged together

- Lower1% = 1th percentile pH value

Above was repeated where SAN had the same value as LOM - since I don't have SAN data I have to approximate it and LOM is the closest by far. I kept freq 0 because there was a separate study that found it never to go that low...

NEW number of outliers:

Still, I got a very few outliers. (max 10). This was all tested with k=1.

## Fix idea number 3:

Looking at the p-value distribution, I realised it looked very conservative, according to [this blogpost](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/) and [this tutorial](https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html) as well as [this tutorial](http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html).

My original p-value distribution:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/LFMM/LFMM_ridge_pvalues_k_1_1.png" width="400" />

“An essential point to understand here is that the FDR is predicated on the “ideal” p-value distribution (flat with a peak at 0)”

This shape of the p-values could be the result of a high GIF (genomic inflation factor). In my case, this was calculated to be k=1 -> 0.6017833, k=2 -> 0.608143, k=7 -> 0.6418904

Manually changing the GIF:

```R
zscore <- pv$score[,1]
new.gif1 <- 0.45
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)
```

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/LFMM/LFMM_ridge_pvalues_k_1_1_adjgif.png" width="400" />

-> Looks much better.

“High genomic inflation factors are caused by population stratification, strong linkage disequilibrium (LD) between SNPs” https://onlinelibrary.wiley.com/doi/full/10.1111/jbg.12419

-> Makes sense that a slightly lower GIF is better for my data. I repeated LFMM tests with the 8 env variables in the table above with lower GIF:

NEW number of outliers (all k=1): (qval < 0.1)

|                 |         | with SAN = LOM | with SAN 8 |
| --------------- | ------- | -------------- | ---------- |
| no date restict | freq7.8 | 869            | 869        |
| no date restict | mins    | 2              | 0          |
| no date restict | ave100  | 1              | 0          |
| no date restict | lower1% | 5              | 0          |
| date restrict   | freq7.8 | 1642           | 1642       |
| date restrict   | mins    | 654            | 10         |
| date restrict   | ave100  | 1110           | 132        |
| date restrict   | lower1% | 768            | 200        |

How does changing K influence the above? (All SAN=LOM, qval<0.1)

|                 |         | k=2  | k=7  |
| --------------- | ------- | ---- | ---- |
| no date restict | freq7.8 | 1015 | 1546 |
| no date restict | mins    | 0    | 35   |
| no date restict | ave100  | 1    | 16   |
| no date restict | lower1% | 5    | 44   |
| date restrict   | freq7.8 | 1635 | 2893 |
| date restrict   | mins    | 607  | 2371 |
| date restrict   | ave100  | 1227 | 2558 |
| date restrict   | lower1% | 620  | 2058 |

Shape of pval distribution for k=7? -> same as k=1

## Fix idea number 4:

Increase minMAF, first from 0.025 to 0.05 -> 588,123 sites, combined with vcf -> 588,119 (extra SNP filtering by angsd)

SAM=LOM

|                 |         | k=1  | k=7  |
| --------------- | ------- | ---- | ---- |
| no date restict | freq7.8 | 647  | 1004 |
| no date restict | mins    | 1    | 3    |
| no date restict | ave100  | 0    | 1    |
| no date restict | lower1% | 2    | 16   |
| date restrict   | freq7.8 | 1090 | 1669 |
| date restrict   | mins    | 445  | 1166 |
| date restrict   | ave100  | 671  | 1263 |
| date restrict   | lower1% | 462  | 1102 |

then to 0.1 -> 337,325 sites, combined with vcf -> 337,325

|                 |         | k=1  | k=7  |
| --------------- | ------- | ---- | ---- |
| no date restict | freq7.8 | 399  | 608  |
| date restrict   | freq7.8 | 646  | 992  |

# Results:

Looking at the 8 environmental variables standardized, such that SD=1, mean=0:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/env_vars_ph.png" width="400" />

As expected mins, bottom 100, and lower 1% are highly correlated with each other. Looking at LFMM outliers (q<0.1, k=7) refects this:

No date restriction:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/image-20220522140014699.png" width="400" />

Date restriction:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/image-20220522140246157.png" width="400" />

With frequency of pH < 7.8:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/image-20220522140646744.png" width="400" />

Date restirced had many more outliers:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/image-20220522140835059.png" width="400" />

And k7 had many more outliers than k1; all k1 outliers also k7 outliers:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/image-20220522140946141.png" width="400" />

In conclusion, choosing k7 and date restricted freq < 7.8 (i.e. qvals_k_7_5_adjgif.csv). 2893 outliers.

```bash
awk -F"," '{print $2}' qvals_k_7_5_adjgif.csv > temp
paste pos_info temp > temp2
awk '$2<0.1' temp2 > temp3
sed -i 's/\.1/\.1_/g' temp3
awk '{print $1}' temp3 > temp4
sed -i 's/^/"/' temp4
sed -i 's/$/"/' temp4
```

Annotated outliers, was at the chi-squared step, and I wanted to integrate the SNP distribution info into it but realised I can only do that if I do same annotation (with enhancers, lncRNA etc) as I did with outliers SNP - previously I only computed annotated01.csv and promoters01.csv so now I am doing this in annotate_outs_5000/annot_outliers/annotate_all.py





TODO: Test temp data too



