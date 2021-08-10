# Loading libraries #
library(vegan)

# Importing the dataset #
taxonomic.div=read.table("results/Taxonomic Diversity.txt",header=T,row.names=1)
taxonomic.comp=read.table("data/ab.rel_orders.txt",header=T,row.names=1)
functional.div=read.table("results/Functional diversity.txt",header=T,row.names=1)
functional.comp=read.table("results/Functional composition.txt",header=T,row.names=1)

# Procrustes Taxonomic and Functional diversity #
pcaTax=prcomp(decostand(taxonomic.div,method="standardize"))
biplot(pcaTax)
pcaFunc=prcomp(decostand(functional.div,method="standardize"))
biplot(pcaFunc)

procrustes(pcaTax,pcaFunc)
protest(pcaTax,pcaFunc)


#dados composição
pcaTax2=prcomp(decostand(taxonomic.comp,method="standardize"))
biplot(pcaTax2)
pcaFunc2=prcomp(decostand(functional.comp,method="standardize"))
biplot(pcaFunc2)

procrustes(pcaTax2,pcaFunc2)
protest(pcaTax2,pcaFunc2)
