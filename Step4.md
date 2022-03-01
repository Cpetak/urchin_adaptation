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
sed 's/0\/0/0/g' cut_filtered_vcf_t | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' | sed 's/NA/9/g' > lfmm_input
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

k=2, binary -> no significant qval

k=2, continous -> 7 under 0.05, 29 under 0.1, for pval 7681 under 0.01, 1242 under 0.001

k=1, continous -> for qval 8 under 0.05, 25 under 0.1, for pval 7490 under 0.01, 1149 under 0.001

so I am choosing continous, and k=1 or k=2 doesn't matter, p-values are very similary, qqplots nearly identical. No k=7 as it doesn't make sence based on PCA plot and also qqplot looks off

K=2 cont it is, p=0.001 is the cut-off

