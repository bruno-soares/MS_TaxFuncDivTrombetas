# Loading libraries #
library(ade4)

# Importing the dataset #
traits=read.table("data/traits.txt",header=T,row.names=1)

# Building a Gower distance matrix #
tabDiet=prep.fuzzy(traits[,1:7],col.blocks=7)
tabMH=prep.fuzzy(traits[,8:12],col.blocks=5)
Ecomor=data.frame(traits[,13:21])
traits.ktab<-ktab.list.df(list(tabDiet,tabMH,Ecomor))
dist.func=dist.ktab(traits.ktab,type=c("F","F","Q"),option=c("scaledBYrange"))
is.euclid(dist.func)

# Measuring the quality of the functional space #
#################################################################################################################################
## R function for computing the quality of functional dendrogramm and multidimensional functional spaces                       ##
##    This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)    ##
#################################################################################################################################

quality_funct_space_fromdist <- function( dist_funct,  nbdim=7,   plot="quality_funct_space") 
{
  
  #loading required libraries
  require(ape)
  require(clue)
  require(cluster)
  require(geometry)
  require(gtools)
  
  ################################################################################################################################# 
  
  # checking data
  if ( ("dist" %in% class(dist_funct) ) ==FALSE )   {  stop(" 'dist_funct' must be of class 'dist'")     }
  if (length(dist_funct)<3)   {  stop(" there must be at least 3 species in 'dist_funct' ")     }
  if (sum(is.na(dist_funct))!=0)   {  stop(" NA are not allowed in 'dist_funct' ")     }
  if (nbdim<2)   {  stop(" 'nbdim' must be higher than 1")     }
  
  # functional distance 
  mat_dissim<-dist_funct
  
  # species names
  nm_sp<-row.names(as.matrix(dist_funct))
  
  ################################################################
  # computing PCoA
  mat_pcoa<-pcoa(mat_dissim)
  
  # changing number of dimensions given number of positive eigenvalues
  nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )
  
  # keeping species coordoinates on the 'nbdim' axes
  mat_coord<-mat_pcoa$vectors[,1:nbdim]
  row.names(mat_coord)<-nm_sp
  colnames(mat_coord)<-paste("PC",1:nbdim,sep="")
  
  # lists to store distance matrices
  dist_raw<-list()
  dist_st<-list()
  
  # computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces 
  for (k in 2:nbdim) 
  {
    eval(parse(text=paste("dist_",k,"D<-dist(mat_coord[,1:",k,"],method='euclidean')", sep="")))
    eval(parse(text=paste("dist_raw$m_",k,"D<-dist_",k,"D", sep="")))
  } # end of k
  
  ################################################################
  # computing an UPGMA-dendrogram then  cophenetic distances between species on this tree
  alg_best_tree<-"UPGMA"
  best_tree<-NA
  dist_best_tree<-NA
  
  best_tree <-hclust( mat_dissim , method="average")
  
  dist_best_tree<-cl_ultrametric(best_tree)
  
  eval(parse(text=paste("dist_raw$t_",alg_best_tree,"<-dist_best_tree", sep="")))
  
  ################################################################ 
  # computing mean squared deviation between initial distance and standardized final distance in the functional space
  meanSD<-rep(NA,nbdim) ; names(meanSD)<-c(paste("t_",alg_best_tree,sep="") ,paste("m_",2:nbdim,"D",sep=""))
  
  x<-mat_dissim # initial distance
  S<-length(nm_sp) # species richness
  
  # for tree
  y<-dist_best_tree
  yst<- y/max(y) * max(x)
  eval(parse(text=paste("dist_st$t_",alg_best_tree,"<-yst", sep="")))
  meanSD[paste("t_",alg_best_tree,sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  
  
  # for muldimensionnal spaces
  for (k in 2:nbdim)  
  {
    eval(parse(text=paste("y<-dist_",k,"D",sep="")))
    yst<- y/max(y) * max(x)
    eval(parse(text=paste("dist_st$m_",k,"D<-dist_",k,"D", sep="")))
    meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k
  
  # list of outputs
  res<-list(meanSD=meanSD, details_funct_space=list( mat_coord=mat_coord, upgma_tree=best_tree, dist_raw=dist_raw, dist_st=dist_st )  )
  
  ################################################################################################################################
  # GRAPHICS if plot has a name
  ################################################################################################################################
  
  if (is.na(plot)==FALSE)
  {
    # window
    if (nbdim<=3) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=600)
      layout(matrix(c(1:(nbdim+1),rep(0,3-nbdim)),1,4,T)) ; layout.show(nbdim+1)  }
    if (nbdim>3 & nbdim<=7) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=1200)
      layout(matrix(c(1:(nbdim+1),rep(0,7-nbdim)),2,4,T)) ; layout.show(nbdim+1) }
    if (nbdim>7 & nbdim<=11) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=1800)
      layout(matrix(c(1:(nbdim+1),rep(0,11-nbdim)),3,4,T)) ; layout.show(nbdim+1) }
    if (nbdim>11 & nbdim<=15) {jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=2400)
      layout(matrix(c(1:(nbdim+1),rep(0,15-nbdim)),4,4,T)) ; layout.show(nbdim+1) }
    if (nbdim>15) { jpeg(paste(plot,".jpeg",sep=""), res=300, width=2400, height=2400)
      layout(matrix(1:16,4,4,T)) ; layout.show(16) ; meanSD_plot<-meanSD[1:15]}	
    
    par(mar=c(4,4,3,3))
    
    # plotting change in meanSD with increasing number of dimensions  
    barplot(height=meanSD,names.arg=names(meanSD), xlab="Functional space", ylab= "Quality (Mean SD)", 
            space=0, cex.names=0.7, col=c("red", rep("blue",nbdim-1) ) )
    
    # plotting quality of each functional space
    
    # functional distances
    x<-as.dist( mat_dissim) 
    
    # dendrogram quality
    eval(parse(text=paste("yst<-dist_st$t_",alg_best_tree, sep="")))
    plot(x,yst , xlab="Initial distance", ylab= "Cophenetic distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
         pch=21, col="red", bg="red", cex=0.3, cex.axis=0.8, cex.lab=0.9 )	
    abline(a=0,b=1)
    title(main=paste(alg_best_tree, "   mSD=",round(meanSD[paste("t_",alg_best_tree,sep="")],4),sep=""), cex.main=1.1, col.main="red", line=0.5   )
    
    # multidimensional spaces quality
    for (k in 2:min(nbdim,15) )  
    {
      eval(parse(text=paste("yst<-dist_st$m_",k,"D",sep="")))
      plot(x,yst, xlab="Initial distance", ylab= "Euclidean distance", xlim=c(0,max(x)), ylim=c(0,max(yst)), 
           pch=21, col="blue", bg="blue", cex=0.3, cex.axis=0.8, cex.lab=0.9  )	
      abline(a=0,b=1)
      title(main=paste(paste(k,"D",sep=""),"   mSD=",round(meanSD[paste("m_",k,"D",sep="")],4),sep=""), cex.main=1.1, col.main="blue", line=0.5  )
    }  # end of k	
    
    
    graphics.off()
    
  } # end of of plot
  
  ################################################################################################################################
  ################################################################################################################################
  
  invisible(res)
  
} # end of function quality_funct_space_fromdist

QFS<-quality_funct_space_fromdist(dist.func, nbdim=10)
QFS$meanSD
which.min(QFS$meanSD)
png("figures/Supplementary Figure 2.png")
plot(1:10,QFS$meanSD[1:10],pch=16,cex=1.3,main="Quality of the functional space",xlab="Axes",ylab="meanSD")
abline(h=0.01,col="gray",lwd=2,lty=2)
dev.off()

coords_pcoa=data.frame(QFS$details_funct_space$mat_coord[,1:3])
envfit_pcoa<-envfit(coords_pcoa,traits,choices=c(1:3))
pcoa(dist.func)

write.table(cbind(envfit_pcoa$vectors$arrows,envfit_pcoa$vectors$r,envfit_pcoa$vectors$pvals),
            "results/Supplementary Table 2.txt")

png("figures/Supplementary Figure 3.png",height=30,width=12,units="cm",res=600)
par(mfrow=c(3,1),mar=c(4,4,2,2))
plot(coords_pcoa$PC2~coords_pcoa$PC1,col="white",xlab="PCoA 1 (38.81%)",ylab="PCoA 2 (15.74%)")
text(coords_pcoa$PC2~coords_pcoa$PC1,labels=rownames(coords_pcoa),cex=1)
plot(coords_pcoa$PC3~coords_pcoa$PC1,col="white",xlab="PCoA 1 (38.81%)",ylab="PCoA 3 (12.31%)")
text(coords_pcoa$PC3~coords_pcoa$PC1,labels=rownames(coords_pcoa),cex=1)
plot(coords_pcoa$PC3~coords_pcoa$PC2,col="white",xlab="PCoA 2 (15.74%)",ylab="PCoA 3 (12.31%)")
text(coords_pcoa$PC3~coords_pcoa$PC2,labels=rownames(coords_pcoa),cex=1)
dev.off()