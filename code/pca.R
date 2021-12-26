C <- as.matrix(read.table("PCangsd_selection.cov"))
ids <- read.table("pca_pops3") #text file with 20 lines of the single word BOD, then 20 lines of CAP etc in the order they appeared in all_rmdups_jo.txt
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
library(tidyverse)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,3], PC2 = e$vectors[,4])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()


for (val in seq(5)){
    for (val2 in seq(5)) {
      df <- data.frame(pop = ids$V1, PC1 = e$vectors[,val], PC2 = e$vectors[,val2])
      df= rownames_to_column(df)
      print(ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
        geom_point(size=3, shape=21) +
        xlab(val) +
        ylab(val2) +
        theme_bw())
    }
}
