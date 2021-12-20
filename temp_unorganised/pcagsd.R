# python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/angsd_noout/allpopstrict_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/angsd_noout/testing_PCangsd_selection -selection -sites_save

# instead of -selection could do -pcadapt

library(RcppCNPy) # Numpy library for R

## function for QQplot
qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

### read in seleciton statistics (chi2 distributed)
s<-npyLoad("testing_PCangsd_selection.selection.npy")

## make QQ plot to QC the test statistics
pdf(file = "QQplot.pdf")
qqchi(s)
dev.off()

# convert test statistic to p-value
pval<-1-pchisq(s,1)

## read positions (hg38)
p<-read.table("testing_PCangsd_selection.sites",colC=c("factor","factor","integer"),sep="_")

names(p)<-c("nothing","chr","pos")

newdf<-cbind(p, pval)

write.csv(newdf,"pcangsd_results.csv")


## make manhatten plot
# plot(-log10(pval),col=p$chr,xlab="Chromosomes",main="Manhattan plot")

## zoom into region
 #w<-range(which(pval<1e-7)) + c(-100,100)
 #keep<-w[1]:w[2]
 #plot(p$pos[keep],-log10(pval[keep]),col=p$chr[keep],xlab="HG38 Position chr2")

## see the position of the most significant SNP
 #p$pos[which.max(s)]
