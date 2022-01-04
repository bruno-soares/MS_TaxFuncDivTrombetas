# Loading libraries #
library(adegenet)
library(MKinfer)

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

# nonparametric bootstrap test for variance between
# pristine and disturbed sites
boot.t.test(abiot$Temperature[1:5],abiot$Temperature[6:9])
boot.t.test(abiot$O2[1:5],abiot$O2[6:9])
boot.t.test(abiot$pH[1:5],abiot$pH[6:9])
boot.t.test(abiot$Conductivity[1:5],abiot$Conductivity[6:9])
boot.t.test(abiot$Turbidity[1:5],abiot$Turbidity[6:9])
boot.t.test(abiot$Av.Widht[1:5],abiot$Av.Widht[6:9])
boot.t.test(abiot$Max.Depth[1:5],abiot$Max.Depth[6:9])
boot.t.test(abiot$Av.Depth[1:5],abiot$Av.Depth[6:9])
boot.t.test(abiot$Streamflow[1:5],abiot$Streamflow[6:9])
boot.t.test(abiot$WaterVelocity[1:5],abiot$WaterVelocity[6:9])
boot.t.test(abiot$CanopyCover[1:5],abiot$CanopyCover[6:9])
boot.t.test(abiot$Subs.SandGravel[1:5],abiot$Subs.SandGravel[6:9])
boot.t.test(abiot$Subs.Leaves[1:5],abiot$Subs.Leaves[6:9])
boot.t.test(abiot$Subs.TrBr[1:5],abiot$Subs.TrBr[6:9])
boot.t.test(abiot$Subs.Vegetation[1:5],abiot$Subs.Vegetation[6:9])
