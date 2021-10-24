# Loading libraries #
library(adegenet)
library(lawstat)

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

# Levene test for heterogeneity of variances
levene.test(abiot$Temperature,groups$Status)
levene.test(abiot$O2,groups$Status)
levene.test(abiot$pH,groups$Status)
levene.test(abiot$Conductivity,groups$Status)
levene.test(abiot$Turbidity,groups$Status)
levene.test(abiot$Av.Widht,groups$Status) #significant
levene.test(abiot$Max.Depth,groups$Status)
levene.test(abiot$Av.Depth,groups$Status)
levene.test(abiot$Streamflow,groups$Status)
levene.test(abiot$WaterVelocity,groups$Status)
levene.test(abiot$CanopyCover,groups$Status)
levene.test(abiot$Subs.SandGravel,groups$Status)
levene.test(abiot$Subs.Leaves,groups$Status)
levene.test(abiot$Subs.TrBr,groups$Status)
levene.test(abiot$Subs.Vegetation,groups$Status)

# t-tests for variance between pristine and disturbed sites
t.test(abiot$Temperature~groups$Status,var.equal=TRUE)
t.test(abiot$O2~groups$Status,var.equal=TRUE)
t.test(abiot$pH~groups$Status,var.equal=TRUE) #marginally significant
t.test(abiot$Conductivity~groups$Status,var.equal=TRUE)
t.test(abiot$Turbidity~groups$Status,var.equal=TRUE) #marginally significant
t.test(abiot$Av.Widht~groups$Status,var.equal=FALSE)
t.test(abiot$Max.Depth~groups$Status,var.equal=TRUE)
t.test(abiot$Av.Depth~groups$Status,var.equal=TRUE)
t.test(abiot$Streamflow~groups$Status,var.equal=TRUE)
t.test(abiot$WaterVelocity~groups$Status,var.equal=TRUE)
t.test(abiot$CanopyCover~groups$Status,var.equal=TRUE)
t.test(abiot$Subs.SandGravel~groups$Status,var.equal=TRUE) #significant
t.test(abiot$Subs.Leaves~groups$Status,var.equal=TRUE) #marginally significant
t.test(abiot$Subs.TrBr~groups$Status,var.equal=TRUE)
t.test(abiot$Subs.Vegetation~groups$Status,var.equal=TRUE)
