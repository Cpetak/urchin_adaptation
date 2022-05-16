#devtools::install_github("bcm-uga/lfmm")
library(lfmm)
library(tseries)

#data("example.data")
#data("skin.exposure")

mydata<-read.matrix("forlapo4", header = FALSE, sep = " ", skip = 0)

#BOD,CAP,FOG,KIB,LOM,SAN,TER
test<-c(2,2,2,0,0,0,0)
test_pos<-rep(test, each=20)
mydata2 <- cbind(mydata, test_pos)

# get genotype data, columns are positions, rows are individuals
#Y <- example.data$genotype
Y <- mydata2

#principal component analysis (PCA) can reveal some ‘structure’ in the genotypic data. 
#We perfom PCA by using the prcomp function as follows.
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
#For the A. thaliana samples, the screeplot indicates that there are around K=6 main components in the data. 
#We will use K=6 latent factors in subsequent analyses of the A. thaliana example.

#The R package contains two main functions for estimating the parameters of LFMMs: 
# ridge_lfmm and lasso_lfmm. 
# The ridge estimates are based on minimimizing a regularized least-squares problem with an L2 penalty.
#Y <- example.data$genotype
#X <- example.data$phenotype #shape: number of individuals x 1

#order of vcf: BOD,CAP,FOG,KIB,LOM,SAN,TER
BOD <- rep(0.05026657, each=20)
CAP <- rep(0.07188418, each=19)
FOG <- rep(0.1780976, each=18)
rest <- rep(c(0.003375068,0.00729927,0,0), each=20)
X<- c(BOD, CAP, FOG, rest)

## Fit an LFMM, i.e, compute B, U, V estimates
# ridge vs lasso:
# ridge minimizes sum of squared residuals + lambda * slope squared -> less sensitive to changes in x, lambda is a scaler, the larger lambda the less sensitive y is to x (i.e. flatter the fitted line)
# overall shrinks parameters, making our predictions less sensitive to them
# vs lasso regression: sum of squared residuals + lambda * absolute_value(slope), again making prediction (y) less sensitive to x
# lasso can shrink parameters to 0, ridge can't (when increasing lambda) -> lasso can exclude useless variables from our model
# ridge is more useful when most parameters are useful -> are most of my parameters useful? I think so!
# so I need ridge... let's do both tho and compare results?
# ridge is also the one used in most papers citing the lfmm paper
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 7) #determined from PCA

#The ridge_lfmm function returns an object that contains the latent variable score matrix U, 
# the latent variable loading matrix U, and B the effect sizes for all SNPs.

## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif") #gif is default, and what others use

#The histogram of test significance values is expected to be flat, with a peak near zero. 
#A QQ-plot is displayed as follows.

pvalues <- pv$calibrated.pvalue 
write.csv(pvalues,"test_pvals")

qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

## Manhattan plot
plot(-log10(pvalues), 
     pch = 19, 
     cex = .9, 
     xlab = "SNP", ylab = "-Log P",
     col = "black")
#points(example.data$causal.set, 
#-log10(pvalues)[example.data$causal.set], 
#type = "h", 
#col = "blue") #The vertical bars indicate the position of the (known) causal loci used to generate the simulated phenotypes.