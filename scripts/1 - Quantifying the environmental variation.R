# Loading libraries #
library(adegenet)

# Importing the dataset #
abiot=read.table("data/abiotic variables.txt",header=T,row.names=1)
groups=read.table("data/grouping variables.txt",header=T,row.names=1) 

#DAPC - Discriminant Analysis of Principal Components (DAPC)
dapc=dapc(decostand(abiot, method="standardize"), groups$Status, center=FALSE,scale=FALSE)
summary(dapc)
dapc$var.contr
write.table(dapc$var.contr,"results/Supplementary Table 1.txt")
dapc$ind.coord

png("figures/Supplementary Figure 1.png")
scatter(dapc,scree.da=FALSE,bg="white",pch=18,solid=.4,cex=3,clab=0,leg=TRUE)
dev.off()
