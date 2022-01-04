# Loading libraries #
library(vegan)
library(SYNCSA)
library(dendextend)
library(lawstat)
library(ggplot2)
library(gridExtra)
library(MKinfer)

# Importing the dataset #
traits=read.table("data/traits.txt",header=T,row.names=1)
abund=read.table("data/cpue_sRec.txt",header=T,row.names=1)
groups=read.table("data/grouping variables.txt",header=T,row.names=1)
abund_ord=read.table("data/ab.rel_orders.txt",header=T,row.names=1)

# taxonomic diversity #
S=specnumber(abund)
H=diversity(abund)
taxonomic=cbind(S,H)
write.table(taxonomic,"results/Taxonomic Diversity.txt")
GerTax2=cbind(taxonomic,groups)

# Nonparametric bootstrap t-test
boot.t.test(GerTax2$S[1:5],GerTax2$S[6:9],var.equal=TRUE)
boot.t.test(GerTax2$H[1:5],GerTax2$H[6:9],var.equal=TRUE)


# Taxonomic composition #
abdOrd=cbind(abund_ord, groups)

# Levene test for taxonomic composition
boot.t.test(abdOrd$X.Cha[1:5],abdOrd$X.Cha[6:9])
boot.t.test(abdOrd$X.Cyp[1:5],abdOrd$X.Cha[6:9])
boot.t.test(abdOrd$X.Per[1:5],abdOrd$X.Cha[6:9])
boot.t.test(abdOrd$X.Sil[1:5],abdOrd$X.Cha[6:9])
boot.t.test(abdOrd$X.Gym.Syn[1:5],abdOrd$X.Cha[6:9])

par(mfrow=c(2,3))
boxplot(X.Cha~factor(Status), data=abdOrd, ylab="%Cha", outline=F, cex.axis=1.5, cex.lab=1.5)
boxplot(X.Cyp~factor(Status), data=abdOrd, ylab="%Cyp", outline=F, cex.axis=1.5, cex.lab=1.5)
boxplot(X.Gym.Syn~factor(Status), data=abdOrd, ylab="%Gym+Syn", outline=F, cex.axis=1.5, cex.lab=1.5)
boxplot(X.Per~factor(Status), data=abdOrd, ylab="%Cic", outline=F, cex.axis=1.5, cex.lab=1.5)
boxplot(X.Sil~factor(Status), data=abdOrd, ylab="%Sil", outline=F, cex.axis=1.5, cex.lab=1.5)
dev.off()



# Functional diversity #
###########################################################################################################################################
# 'multidimFD': function to compute and illustrate multidimensional functional diversity indices for a set of species assemblages
# For details about indices formulas see Mouillot et al. 2013, Trends in ecology and Evolution (28:167-177) and references therein.
# Provided by Sebastien Villéger at http://villeger.sebastien.free.fr/                                                                  
##########################################################################################################################################

multidimFD<-function(coord, weight, check_species_pool=TRUE, verb=TRUE,
                     folder_plot=NULL, nm_asb_plot=NULL, Faxes_plot=NULL, Faxes_nm_plot=NULL, 
                     plot_pool=TRUE, col_bg="grey90", col_sp_pool="grey30", pch_sp_pool="+", cex_sp_pool=1,
                     pch_sp=21, col_sp="#1145F0", transp=50 ) 
  
{
  
  # library required for indices computation
  require (geometry)
  require(ape)
  
  # saving name of current working directory
  current_wd<-getwd()
  
  ##############################################################################
  # checking inputs
  
  
  # coordinates of species in the functional space
  if( nrow(coord)<2 ) stop(paste(" error: there must be at least 2 species in the dataset"))
  if( ncol(coord)<2 ) stop(paste(" error: there must be at least 2 functional axes"))
  if( is.numeric(coord)==FALSE ) stop(paste(" error: 'coord' is not numeric"))
  if( is.na(sum(coord)) ) stop(paste(" error: NA in 'coord'"))
  if ( is.null(colnames(coord)) ) { colnames(coord)<-paste("Axis", 1:ncol(coord),sep="_") } # if no column names in 'coord' default value
  
  # dominance of species in assemblages
  if( is.matrix(weight)==FALSE ) stop( " 'weight' should be an object of type matrix")
  if( is.numeric(weight)==FALSE ) stop(paste(" error: 'weight' is not numeric"))
  if( is.na(sum(weight)) ) stop(paste(" error: NA in 'weight'"))
  if( min(weight)<0 ) stop(paste(" error: negative value in 'weight'"))
  if(min(apply(weight,1,sum))==0 ) 
    stop(paste(" error: all rows of 'weight' should have a sum stricly positive, i.e. all assemblage must host at least one species"))
  
  # match between datasets
  if( sum(colnames(weight) == row.names(coord))!= nrow(coord) ) stop(paste(" error: 'weight' does not have the same column names than row names of 'coord'"))
  
  # checking graphical parameters
  if( length(pch_sp)!=1 ) stop(paste(" error:'pch_sp' should contain only one value"))
  if( length(col_sp)!=1 ) stop(paste(" error:'col_sp' should contain only one value"))
  if( length(col_bg)!=1 ) stop(paste(" error:'col_bg' should contain only one value"))
  if( length(col_sp_pool)!=1 ) stop(paste(" error:'col_sp_pool' should contain only one value"))
  if( length(pch_sp_pool)!=1 ) stop(paste(" error:'pch_sp_pool' should contain only one value"))
  if( length(cex_sp_pool)!=1 ) stop(paste(" error:'cex_sp_pool' should contain only one value"))
  
  
  # checking species pool
  if (check_species_pool==TRUE)
  {
    if(min(apply(weight,2,sum))==0 ) 
      stop(paste(" error: all columns of 'weight' should have a sum stricly positive, i.e. all species must occur in at least one assemblage"))
  }# end of check species pool
  
  
  ##############################################################################
  # info on study case
  
  # number and names of axes
  nm_axes<-colnames(coord)
  nb_axes<-length(nm_axes)
  
  # number and names of assemblages
  nm_asb<-row.names(weight)
  nb_asb<-length(nm_asb)
  
  # matrix to store results
  indices<-c( "Nb_sp", "Tot_weight", paste("min",nm_axes,sep="_"), paste("max",nm_axes,sep="_"), paste("range",nm_axes,sep="_"), 
              paste("FIde",nm_axes,sep="_"), c("FRic","FDiv","FEve","FDis","FSpe", "FOri") )
  FD<-matrix(NA, nb_asb, length(indices), dimnames=list(nm_asb,indices))
  
  ##############################################################################
  # preliminary computation at the species pool level
  
  #######################################
  # originality of each species: distance to nearest neighbour among the global pool of species
  dist_sp<-as.matrix(dist(coord,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
  orig_sp<-apply(dist_sp, 1, min, na.rm=T )
  # identity of Nearest Neighbour
  NN<-dist_sp ; NN<-NN-apply(NN,1,min,na.rm=T) ; NN[which(NN!=0)]<-NA   ; NN[which(NN==0)]<-1
  
  # specialization of each species: distance to centroid of the global pool of species
  centroid_sp<-apply(coord,2,mean) # coordinates of the center of gravity of the vertices (B)
  spec_sp<-apply(coord, 1, function(x) { (sum((x-centroid_sp)^2) )^0.5} )
  
  # convex hull volume of the species pool
  FRic_pool<-convhulln(coord,"FA")$vol
  
  #######################################
  # setting same graphical parameters for all assemblages
  
  # setting working directory to store jpeg files
  if (is.null(folder_plot) ) {  folder_plot<-current_wd }
  
  # setting folder for saving jpeg files
  test_folder_plot<-try( setwd(folder_plot) , silent=TRUE)
  if ( class(test_folder_plot) =="try-error") {
    folder_plot<-current_wd 
    print(paste(" /!\    WARNING: '",folder_plot," does not exist', jpeg files have been saved in '",current_wd, "'",sep="" ))
  } # end of if
  setwd(folder_plot)
  
  # range of 'coord' extended by 5%
  rge_coord<-range(coord)
  extrge_coord<-c( rge_coord[1]-0.05*(rge_coord[2]-rge_coord[1]) , rge_coord[2]+0.05*(rge_coord[2]-rge_coord[1]))
  
  # pretty axes labels within observed range
  lab<-pretty(extrge_coord, n=4) # labels
  lab<-lab[which(lab>=rge_coord[1] & lab<=rge_coord[2]) ]# filtering
  
  # axes limits a,d range
  Faxes_lim<-extrge_coord
  rge_Faxes_lim<-Faxes_lim[2]-Faxes_lim[1]
  
  ###########################################
  # function to draw functional space given axes limits and background color, option: lengend for species weights
  functional_space<-function(axes_xy, nm_axes_xy, col_bg, plot_pool=plot_pool, legend_weight=FALSE) 
  {
    # setting margins size of plot
    par( pty="s", mar=c(3,3,3,3) ) 
    
    plot(Faxes_lim,Faxes_lim,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=Faxes_lim,ylim=Faxes_lim) # default window
    rect(Faxes_lim[1],Faxes_lim[1],Faxes_lim[2],Faxes_lim[2], col="white")   # customized window
    
    # customized X and Y axes
    for (k in 1:2)
    {
      axis(side=k, at=lab, labels=F, tcl=-0.3, pos=Faxes_lim[1])  # ticks
      mtext(side=k, lab, at=lab, line=-0.2, cex=0.9, las=1) # labels
      mtext(side=k, nm_axes_xy[k], cex=1,line=1.5, font=2) # title  
    } # end of k
    
    
    if (plot_pool==TRUE)
    {
      # grey background
      rect(Faxes_lim[1],Faxes_lim[1],Faxes_lim[2],Faxes_lim[2], col=col_bg)   
      
      # convex hull of species pool
      vert0<-convhulln( coord[,axes_xy] ,"Fx TO 'vert.txt'")
      vert1<-scan("vert.txt",quiet=T)
      vert_ij<-(vert1+1)[-1]
      polygon(coord[vert_ij,axes_xy], border=NA,col="white")   
      
      # all species
      points( coord[,axes_xy[1] ], coord[,axes_xy[2] ] , pch=pch_sp_pool, col=col_sp_pool, cex=cex_sp_pool)
      
    }# end of if plot_pool
    
    
    # legend for abundance 
    if (legend_weight==TRUE)  
    {
      rect(max(Faxes_lim)-0.25*rge_Faxes_lim, min(Faxes_lim), max(Faxes_lim), min(Faxes_lim)+0.12*rge_Faxes_lim, col="white")
      symbols(max(Faxes_lim)-0.19*rge_Faxes_lim, min(Faxes_lim)+0.06*rge_Faxes_lim, circles=sqrt(0.1)*0.075*rge_Faxes_lim, 
              inches=FALSE, bg="black", fg="black", add=TRUE, lwd=1.5)
      text(max(Faxes_lim)-0.15*rge_Faxes_lim, min(Faxes_lim)+0.06*rge_Faxes_lim,"10%", adj=c(0,0.5) ) 
    }# end of if legend
    
  }# end of function functional_space
  ###########################################
  
  # end of preliminary computation
  
  ##############################################################################
  
  # loop on assemblages for computing and plotting functional diversity indices
  for (k in nm_asb)
  {
    
    ###########################################################
    # preparing data
    
    # names, number, weight and coordinates of of species present
    weight_k<-weight[k,]
    nm_sp_k<-row.names(coord)[which(weight_k>0)]
    nb_sp_k<-length(nm_sp_k)
    weight_sp_k<-weight[k,nm_sp_k]
    coord_sp_k<-coord[nm_sp_k,]
    if(nb_sp_k==1) { coord_sp_k<-matrix(coord_sp_k,nrow=1,dimnames=list(nm_sp_k,colnames(coord)) ) } # matrix object
    
    # names of species absent
    nm_sp_absent_k<-names(which(weight[k,]==0))
    
    # total weight
    FD[k,"Tot_weight"]<-sum(weight_sp_k)
    
    #relative weight 
    rel_weight_sp_k<-weight_sp_k/sum(weight_sp_k)
    
    # species richness
    FD[k,"Nb_sp"]<-nb_sp_k
    
    ###########################################################
    # computing indices values on each axis
    
    # range of values
    for (z in nm_axes) 
    {
      FD[k,paste("min",z,sep="_")]<-min(coord_sp_k[,z])
      FD[k,paste("max",z,sep="_")]<-max(coord_sp_k[,z])
      FD[k,paste("range",z,sep="_")]<-FD[k,paste("max",z,sep="_")]-FD[k,paste("min",z,sep="_")]
    }# end of z
    
    # abundance-weighted mean values
    FD[k,paste("FIde",nm_axes,sep="_")]<-rel_weight_sp_k%*%coord_sp_k
    
    ###########################################################  
    # multivariate indices
    
    
    # indices based on vertices and volume of convex hull, only if more species than number of axes
    
    if (nb_sp_k>nb_axes) {
      
      ########################
      # Functional richness = convex hull volume
      FD[k,"FRic"]<-round(convhulln(coord_sp_k,"FA")$vol/FRic_pool,6)
      
      
      ########################
      # Functional Divergence
      
      # identity of vertices
      vert0<-convhulln(coord_sp_k,"Fx TO 'vert.txt'")
      vert1<-scan("vert.txt",quiet=T)
      vertices_k<-(vert1+1)[-1]
      
      # coordinates of the center of gravity of the vertices (B)
      B_k<-apply(coord_sp_k[vertices_k,],2,mean)
      
      # Euclidean dstance to B (dB)
      dB_k<-apply(coord_sp_k, 1, function(x) { (sum((x-B_k)^2) )^0.5} )
      
      # mean of dB values and deviations to mean 
      meandB_k<-mean(dB_k)
      devdB_k<-dB_k-meandB_k
      
      # abundance-weighted mean deviation
      abdev_k<- rel_weight_sp_k*devdB_k
      ababsdev_k<- rel_weight_sp_k*abs(devdB_k)
      
      FD[k,"FDiv"]<-round( (sum(abdev_k)+meandB_k) / (sum(ababsdev_k)+meandB_k) ,6)
          }
    
    ##########################
    # Functional Evenness
    
    if (nb_sp_k>=3) {
      
      # inter-species Euclidean distance
      distT_k<-dist(coord_sp_k, method="euclidian")
      
      # topology of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
      linkmst_k<-mst(distT_k)
      mstvect_k<-as.dist(linkmst_k)
      
      # pairwise cumulative relative abundances and conversion into 'dist' class
      ab2_k<-matrix(0,nrow=nb_sp_k,ncol=nb_sp_k)
      for (q in 1:nb_sp_k)
        for (r in 1:nb_sp_k)
          ab2_k[q,r]<-rel_weight_sp_k[q]+rel_weight_sp_k[r] # end of q,r
      ab2vect_k<-as.dist(ab2_k)
      
      # EW index for the (S-1) segments
      EW_k<-rep(0,nb_sp_k-1)
      flag<-1
      for (m in 1:((nb_sp_k-1)*nb_sp_k/2))
      {if (mstvect_k[m]!=0) {EW_k[flag]<-distT_k[m]/(ab2vect_k[m]) ; flag<-flag+1}}  # end of m
      
      # PEW index and comparison with 1/S-1
      minPEW_k<-rep(0,nb_sp_k-1)  ;  OdSmO_k<-1/(nb_sp_k-1)
      for (l in 1:(nb_sp_k-1))
        minPEW_k[l]<-min( (EW_k[l]/sum(EW_k)) , OdSmO_k )  # end of l
      
      # FEve
      FD[k,"FEve"]<-round( ( (sum(minPEW_k))- OdSmO_k) / (1-OdSmO_k ) ,6)
      
    }
    
    ##########################
    # Functional Dispersion: abundance-weighted mean distance to abundance-weighted centroid 
    # scaled by maximum value possible given species pool (i.e. the two most distant species have half of total weight)
    dist_centr_sp_k<-apply(coord_sp_k, 1, function(x) { (sum((x-FD[k,paste("FIde",nm_axes,sep="_")])^2) )^0.5} ) # distance to abundance-weighted centroid 
    FD[k,"FDis"]<-(rel_weight_sp_k %*% dist_centr_sp_k) / ( max(dist_sp, na.rm=T) /2 )
    
    ##########################
    # functional specialization : abundance-weighted mean distance to centroid of the global pool of species
    # scaled by maximum value possible given species pool (i.e. an assmeblage hosting only the most specialized species)
    FD[k,"FSpe"]<-(rel_weight_sp_k %*% spec_sp[nm_sp_k])/ max(spec_sp)
    
    ##########################
    # functional originality : abundance-weighted mean distance to nearest neighbour in the global pool of species 
    # scaled by maximum value possible given species pool (i.e. an assmeblage hosting only the most original species)
    FD[k,"FOri"]<-(rel_weight_sp_k %*% orig_sp[nm_sp_k])/max(orig_sp)
    
    
    ######################################################################################################################
    # End of indices computation
    ######################################################################################################################
    
    
    ###########################################################
    # if graphical output
    
    if (k %in% nm_asb_plot)
    {
      
      # axes to plot if not specified: up to first 4 axes 
      if ( is.null(Faxes_plot) ) { Faxes_plot<-colnames(coord)[ 1:(min(c(ncol(coord),4)))]  }
      
      # igf not specififed default axes names
      if ( is.null(Faxes_nm_plot) ) { Faxes_nm_plot<-Faxes_plot }
      
      # checking inputs
      if( sum( nm_asb_plot %in% row.names(weight) ) != length(nm_asb_plot) ) stop(paste(" error: 'nm_asb_plot' should be subset of 'coord' row names"))
      
      if( length(Faxes_plot) <2 | length(Faxes_plot) >4) stop(paste("length of 'Faxes_plot' should be 2, 3 or 4 "))
      if( sum( Faxes_plot %in% colnames(coord) ) != length(Faxes_plot) ) stop(paste(" error: 'Faxes_plot' should be subset of 'coord' column names"))
      if( length(Faxes_plot) != length(Faxes_nm_plot) ) stop(paste("length of 'Faxes_plot' should match length of 'Faxes_nm_plot' "))
      
      # shortening object names
      Faxes<-Faxes_plot
      Faxes_nm<-Faxes_nm_plot
      
      # number of axes
      nb_Faxes<-length(Faxes_plot)
      
      ################################
      # loop on pairs of axes
      for (i in 1:(nb_Faxes-1) )
        for (j in (i+1):nb_Faxes )
        {
          
          # creating jpeg file with 8 panels
          nmjpeg<-paste(k,"_",Faxes[i] ,"_",Faxes[j],".jpeg",sep="")
          jpeg(file=nmjpeg, res=150, width=1800, height=1200)
          layout(matrix(c(1:6),2,3,T)) ; layout.show(6)
          
          # setting margins size of plot
          par( pty="s", mar=c(3,3,3,3) ) 
          
          ################################
          # First panel = Functional Identity and Functional dispersion
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          # mean position on each axis 
          segments(FD[k,paste("FIde_",Faxes[i],sep="")],FD[k,paste("FIde_",Faxes[j],sep="")], FD[k,paste("FIde_",Faxes[i],sep="")], min(Faxes_lim), lwd=1.5, col=col_sp, lty=2)
          segments(FD[k,paste("FIde_",Faxes[i],sep="")],FD[k,paste("FIde_",Faxes[j],sep="")], min(Faxes_lim) ,FD[k,paste("FIde_",Faxes[j],sep="")], lwd=1.5, col=col_sp, lty=2)
          points( FD[k,paste("FIde_",Faxes[i],sep="")],FD[k,paste("FIde_",Faxes[j],sep="")], pch=22, bg=col_sp, col=col_sp,cex=2.5)
          
          # distance to abundance weighted centroid
          segments( FD[k,paste("FIde_",Faxes[i],sep="")],FD[k,paste("FIde_",Faxes[j],sep="")], coord_sp_k[,Faxes[i]], coord_sp_k[,Faxes[j]],col=col_sp, lty=1, lwd=2)
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab_k<-sqrt(rel_weight_sp_k)*0.075*rge_Faxes_lim
          o_weight<-order(rel_weight_sp_k, decreasing = TRUE)
          symbols(coord_sp_k[o_weight,Faxes[i]], coord_sp_k[o_weight,Faxes[j]], circles=sizeab_k[o_weight], inches=FALSE, 
                  bg=col_sp, fg="black", add=TRUE)
          
          # FDis and FIde values
          mtext(side=3, paste("FDis=",round(FD[k,'FDis'],3),sep=""), at=mean(Faxes_lim), line=0.8, cex=1.1,adj=0.5, font=2)               
          mtext(side=3, paste("FIde(x)=",round(FD[k,paste("FIde_",Faxes[i],sep="")],3),sep=""), at=min(Faxes_lim), line=-0.4, cex=0.9,adj=0, font=2)     
          mtext(side=3, paste("FIde(y)=",round(FD[k,paste("FIde_",Faxes[j],sep="")],3),sep=""), at=max(Faxes_lim), line=-0.4, cex=0.9,adj=1, font=2)     
          
          ################################
          # Second panel = functional richness
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          # range on each axis
          dec1<-rge_Faxes_lim *0.02
          segments( FD[k,paste("min_",Faxes[i],sep="")], min(Faxes_lim)-dec1, FD[k,paste("min_",Faxes[i],sep="")], min(Faxes_lim)+dec1, col=col_sp , lwd=3) # min x
          segments( FD[k,paste("max_",Faxes[i],sep="")], min(Faxes_lim)-dec1, FD[k,paste("max_",Faxes[i],sep="")], min(Faxes_lim)+dec1, col=col_sp , lwd=3) # max x
          segments( min(Faxes_lim)-dec1, FD[k,paste("min_",Faxes[j],sep="")], min(Faxes_lim)+dec1, FD[k,paste("min_",Faxes[j],sep="")],  col=col_sp , lwd=3) # min y
          segments( min(Faxes_lim)-dec1, FD[k,paste("max_",Faxes[j],sep="")], min(Faxes_lim)+dec1, FD[k,paste("max_",Faxes[j],sep="")],  col=col_sp , lwd=3) # max y
          
          
          # projected convex hull in 2D
          if (nb_sp_k>nb_axes) {
            
            vert0<-convhulln(coord_sp_k[,Faxes[c(i,j)]],"Fx TO 'vert.txt'")
            vert1<-scan("vert.txt",quiet=T) ; vertices2D<-(vert1+1)[-1]
            polygon(coord_sp_k[vertices2D,Faxes[c(i,j)]],border=NA,col=paste(col_sp,transp,sep=""))
            
            # all points (empty) then filling points being vertices in nD
            points(coord_sp_k[,Faxes[i]], coord_sp_k[,Faxes[j]], pch=21, bg="white", col="black",cex=2)
            points(coord_sp_k[vertices_k,Faxes[i]], coord_sp_k[vertices_k,Faxes[j]], pch=21,bg=col_sp, col="black",cex=2)
            
            # index    
            mtext(side=3, paste("FRic=",round(FD[k,'FRic'],3),sep=""), at=mean(Faxes_lim), line=0.8, cex=1.1,adj=0.5, font=2)  
          }# end of FRic computable
          
          mtext(side=3, paste(Faxes_nm[i], " [",round(FD[k,paste("min_",Faxes[i],sep="")],1),";",round(FD[k,paste("max_",Faxes[i],sep="")],1),"]",sep=""),
                at=min(Faxes_lim), line=-0.4, cex=0.8,adj=0) 
          mtext(side=3, paste(Faxes_nm[j], " [",round(FD[k,paste("min_",Faxes[j],sep="")],1),";",round(FD[k,paste("max_",Faxes[j],sep="")],1),"]",sep=""),
                at=max(Faxes_lim), line=-0.4, cex=0.8,adj=1) 
          
          ###############################################################
          # Third panel = functional Divergence
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          if (nb_sp_k>nb_axes) {
            # projected convex hull in 2D
            vert0<-convhulln(coord_sp_k[,Faxes[c(i,j)]],"Fx TO 'vert.txt'")
            vert1<-scan("vert.txt",quiet=T) ; vertices2D<-(vert1+1)[-1]
            polygon(coord_sp_k[vertices2D,Faxes[c(i,j)]],border=NA,col=paste(col_sp,transp,sep=""))
            
            # distance to centroid of vertices
            segments( B_k[Faxes[i]], B_k[Faxes[j]], coord_sp_k[,Faxes[i]], coord_sp_k[,Faxes[j]],col=col_sp, lty=1, lwd=2)
            points( B_k[Faxes[i]], B_k[Faxes[j]], pch=23,col=col_sp, bg=col_sp,cex=2.5)
          }# end of FRic computable
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab_k<-sqrt(rel_weight_sp_k)*0.075*rge_Faxes_lim
          o_weight<-order(rel_weight_sp_k, decreasing = TRUE)
          symbols(coord_sp_k[o_weight,Faxes[i]], coord_sp_k[o_weight,Faxes[j]], circles=sizeab_k[o_weight], inches=FALSE, 
                  bg=col_sp, fg="black", add=TRUE)
          
          # FDiv index    
          mtext(side=3, paste("FDiv=",round(FD[k,'FDiv'],3),sep=""), at=mean(Faxes_lim), line=-0.1, cex=1.1,adj=0.5, font=2)               
          
          ###############################################################
          # Fourth panel = functional evenness
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          # MST
          if (nb_sp_k>=3) {
            for (x in 1:nrow(linkmst_k))
              for (y in 1:nrow(linkmst_k))
                if (linkmst_k[y,x]==1 & y>x) {
                  segments(coord_sp_k[y,Faxes[i]], coord_sp_k[y,Faxes[j]],coord_sp_k[x,Faxes[i]], coord_sp_k[x,Faxes[j]], 
                           col=col_sp, lwd=1.5) }# end of if link on MST
          }# end of FEve computable
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab_k<-sqrt(rel_weight_sp_k)*0.075*rge_Faxes_lim
          o_weight<-order(rel_weight_sp_k, decreasing = TRUE)
          symbols(coord_sp_k[o_weight,Faxes[i]], coord_sp_k[o_weight,Faxes[j]], circles=sizeab_k[o_weight], inches=FALSE, 
                  bg=col_sp, fg="black", add=TRUE)
          
          # FEve index    
          mtext(side=3, paste("FEve=",round(FD[k,'FEve'],3),sep=""), at=mean(Faxes_lim), line=-0.1, cex=1.1,adj=0.5, font=2)               
          
          
          ###############################################################
          # Fifth panel = functional specialization
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          # distance to centroid of all points
          segments( centroid_sp[Faxes[i]], centroid_sp[Faxes[j]], coord_sp_k[,Faxes[i]], coord_sp_k[,Faxes[j]],col=col_sp, lty=3, lwd=2)
          points( centroid_sp[Faxes[i]], centroid_sp[Faxes[j]], pch=23,col="black",bg="black",cex=2.5)
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab_k<-sqrt(rel_weight_sp_k)*0.075*rge_Faxes_lim
          o_weight<-order(rel_weight_sp_k, decreasing = TRUE)
          symbols(coord_sp_k[o_weight,Faxes[i]], coord_sp_k[o_weight,Faxes[j]], circles=sizeab_k[o_weight], inches=FALSE, 
                  bg=col_sp, fg="black", add=TRUE)
          
          # FSpe index    
          mtext(side=3, paste("FSpe=",round(FD[k,'FSpe'],3),sep=""), at=mean(Faxes_lim), line=-0.1, cex=1.1,adj=0.5, font=2)               
          
          ###############################################################
          # Sixth panel = functional originality
          functional_space( c(Faxes[i] ,Faxes[j]), c(Faxes_nm[i] ,Faxes_nm[j]), col_bg=col_bg, plot_pool=plot_pool, legend_weight=TRUE )
          
          for (z in row.names(coord_sp_k) )
          {
            nm_NN_z<-names(which(NN[z,]==1)[1])
            arrows( coord_sp_k[z,Faxes[i]], coord_sp_k[z,Faxes[j]], coord[nm_NN_z,Faxes[i]], coord[nm_NN_z,Faxes[j]],
                    col="black", lwd=1.8, length=0.1, angle=20)
          } # end of k
          
          # abundances, scaling: point area proportional to relative abundance, if relab=100%, circle diameter=15% of axis range
          sizeab_k<-sqrt(rel_weight_sp_k)*0.075*rge_Faxes_lim
          o_weight<-order(rel_weight_sp_k, decreasing = TRUE)
          symbols(coord_sp_k[o_weight,Faxes[i]], coord_sp_k[o_weight,Faxes[j]], circles=sizeab_k[o_weight], inches=FALSE, 
                  bg=col_sp, fg="black", add=TRUE)
          
          # FSpe index    
          mtext(side=3, paste("FOri=",round(FD[k,'FOri'],3),sep=""), at=mean(Faxes_lim), line=-0.1, cex=1.1,adj=0.5, font=2)               
          
        }# end of i,j (pair of axes)
      
      # closing jpeg
      graphics.off()
      ################################
      
      
    }# end of plot of FD indices
    ###########################################################
    
    # printing step achieved
    if (verb==TRUE) print(paste("FD of assemblage '",k,"' computed",sep="") )
  }# end of working on assemblage k
  ###########################################################  
  
  # returning to current working directory
  setwd(current_wd)
  
  # returning results	
  return(FD)	
  
}# end of function multidimFD

funcdiv<-multidimFD(QFS$details_funct_space$mat_coord[,1:3],as.matrix(abund))
funcdiv<-funcdiv[,c(15:20)]
funcdiv<-as.data.frame(cbind(funcdiv,FRedundancy))
write.table(funcdiv,"results/Functional diversity.txt")

# Levene test for functional structure
boot.t.test(funcdiv$FRic[1:5],funcdiv$FRic[6:9])
boot.t.test(funcdiv$FDiv[1:5],funcdiv$FDiv[6:9])
boot.t.test(funcdiv$FEve[1:5],funcdiv$FEve[6:9])
boot.t.test(funcdiv$FDis[1:5],funcdiv$FDis[6:9])
boot.t.test(funcdiv$FOri[1:5],funcdiv$FOri[6:9])
boot.t.test(funcdiv$FSpe[1:5],funcdiv$FSpe[6:9])


# Functional composition #
funccomp<-multidimFD(QFS$details_funct_space$mat_coord[,1:3],as.matrix(abund))
funccomp<-as.data.frame(funccomp[,c(12:14)])
write.table(funccomp,"results/Functional composition.txt")

par(mfrow=c(1,3))
boxplot(funccomp$FIde_PC1~groups$Status,ylab="FIde_PC1",cex.axis=2, outline=F)
boxplot(funccomp$FIde_PC2~groups$Status,ylab="FIde_PC2",cex.axis=2, outline=F)
boxplot(funccomp$FIde_PC3~groups$Status,ylab="FIde_PC3",cex.axis=2, outline=F)
dev.off()

# Levene test for functional composition
boot.t.test(funccomp$FIde_PC1[1:5],funccomp$FIde_PC1[6:9])
boot.t.test(funccomp$FIde_PC2[1:5],funccomp$FIde_PC2[6:9])
boot.t.test(funccomp$FIde_PC3[1:5],funccomp$FIde_PC3[6:9])


# Supplementary Figure 4 #
Richness<-ggplot(data=GerTax2)+
  geom_boxplot(mapping=aes(y=S,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=S,x=factor(Status)))+
  ylab("Richness")+
  xlab(NULL)+
  geom_text(aes(y=11.9,x=0.5,label=paste('A'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

Shannon<-ggplot(data=GerTax2)+
  geom_boxplot(mapping=aes(y=H,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=H,x=factor(Status)))+
  ylab("Shannon's diversity")+
  xlab(NULL)+
  geom_text(aes(y=1.72,x=0.5,label=paste('B'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

png("figures/Supplementary Figure 4.png",units="cm",width = 16,height=10,res=600)
grid.arrange(Richness,Shannon,nrow=1)
dev.off()




# Supplementary Figure 5 #
Characiformes<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Cha,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Cha,x=factor(Status)))+
  ylab("% Characiformes")+
  xlab(NULL)+
  geom_text(aes(y=52,x=0.52,label=paste('A'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

Cyprinodontiformes<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Cyp,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Cyp,x=factor(Status)))+
  ylab("% Cyprinodontiformes")+
  xlab(NULL)+
  geom_text(aes(y=69,x=0.52,label=paste('B'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

Gymnotiformes<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Gym.Syn,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Gym.Syn,x=factor(Status)))+
  ylab("% Gymnotiformes & Synbranchiformes")+
  xlab(NULL)+
  geom_text(aes(y=9,x=0.52,label=paste('C'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

Cichliformes<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Per,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Per,x=factor(Status)))+
  ylab("% Cichliformes")+
  xlab(NULL)+
  geom_text(aes(y=42,x=0.52,label=paste('D'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

Siluriformes<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Sil,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Sil,x=factor(Status)))+
  ylab("% Siluriformes")+
  xlab(NULL)+
  geom_text(aes(y=17,x=0.52,label=paste('E'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

png("figures/Supplementary Figure 5.png",units="cm",width=16,height=12,res=600)
grid.arrange(Characiformes,Cyprinodontiformes,Gymnotiformes,Cichliformes,
             Siluriformes)
dev.off()


# Supplementary Figure 6 #
fric<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FRic,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FRic,x=factor(groups$Status)))+
  ylab("Functional Richness")+
  xlab(NULL)+
  geom_text(aes(y=0.56,x=0.52,label=paste('A'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

fdiv<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FDiv,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FDiv,x=factor(groups$Status)))+
  ylab("Functional Divergence")+
  xlab(NULL)+
  geom_text(aes(y=0.874,x=0.52,label=paste('B'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

feve<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FEve,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FEve,x=factor(groups$Status)))+
  ylab("Functional Evenness")+
  xlab(NULL)+
  geom_text(aes(y=0.75,x=0.52,label=paste('C'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

fdis<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FDis,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FDis,x=factor(groups$Status)))+
  ylab("Functional Dispersion")+
  xlab(NULL)+
  geom_text(aes(y=0.6,x=0.52,label=paste('D'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

fori<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FOri,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FOri,x=factor(groups$Status)))+
  ylab("Functional Originality")+
  xlab(NULL)+
  geom_text(aes(y=0.357,x=0.52,label=paste('E'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

fspe<-ggplot()+
  geom_boxplot(mapping=aes(y=funcdiv$FSpe,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funcdiv$FSpe,x=factor(groups$Status)))+
  ylab("Functional Specialization")+
  xlab(NULL)+
  geom_text(aes(y=0.55,x=0.52,label=paste('F'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

png("figures/Supplementary Figure 6.png",units="cm",width=16,height=18,res=600)
grid.arrange(fric,fdiv,feve,fdis,fori,fspe)
dev.off()



# Supplementary Figure 7 #
PC1<-ggplot()+
  geom_boxplot(mapping=aes(y=funccomp$FIde_PC1,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funccomp$FIde_PC1,x=factor(groups$Status)))+
  ylab("CWM (PC1)")+
  xlab(NULL)+
  geom_text(aes(y=0.073,x=0.52,label=paste('A'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

PC2<-ggplot()+
  geom_boxplot(mapping=aes(y=funccomp$FIde_PC2,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funccomp$FIde_PC2,x=factor(groups$Status)))+
  ylab("CWM (PC2)")+
  xlab(NULL)+
  geom_text(aes(y=0.1815,x=0.52,label=paste('B'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

PC3<-ggplot()+
  geom_boxplot(mapping=aes(y=funccomp$FIde_PC3,x=factor(groups$Status)),fill="grey")+
  geom_jitter(mapping=aes(y=funccomp$FIde_PC3,x=factor(groups$Status)))+
  ylab("CWM (PC3)")+
  xlab(NULL)+
  geom_text(aes(y=-0.007,x=0.52,label=paste('C'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 6, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

png("figures/Supplementary Figure 7.png",units="cm",width = 16,height=8,res=600)
grid.arrange(PC1,PC2,PC3,nrow=1)
dev.off()




# Richness effects on FRic and FOri #
png("figures/Supplementary Figure 8.png",units="cm",width = 16,height=10,res=600)
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(funcdiv$FOri~S,xlab="Species Richness",ylab="Functional Originality",
     pch=16,col=alpha(c("red","red","red","red","red","blue","blue","blue",
                        "blue","blue"),0.6))
text(6.2,0.36,paste('B'))
plot(funcdiv$FRic~S,xlab="Species Richness",ylab="Functional Richness",
     pch=16,col=alpha(c("red","red","red","red","red","blue","blue","blue",
                  "blue","blue"),0.6))
text(6.2,0.57,paste('A'))
dev.off()

summary(lm(funcdiv$FRic~S))
summary(lm(funcdiv$FOri~S))


# Figure 2 #
Siluriformes_main<-ggplot(data=abdOrd)+
  geom_boxplot(mapping=aes(y=X.Sil,x=factor(Status)),fill="grey")+
  geom_jitter(mapping=aes(y=X.Sil,x=factor(Status)))+
  ylab("% Siluriformes")+
  xlab(NULL)+
  geom_text(aes(y=17,x=0.52,label=paste('B'),hjust=0.5,vjust=2.1),size=3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # opcoes graficas
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(face = "bold", colour = "black", size = 5))

png("figures/Figure 2.png",units="cm",width = 18,height=10,res=600)
grid.arrange(Richness,Siluriformes_main,nrow=1)
dev.off()
