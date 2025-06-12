
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(lattice)
library(tidyr)
library(SeqArray)
library(stringr)
library(doMC)
registerDoMC(20)
library(lubridate)
library(ggsci)
library(ggpubfigs)
library(adegenet)


genofile <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz.gds" , allow.fork=T)
samps <- read.gdsn(index.gdsn(genofile, "sample.id")) 

#get metadata
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
metadata[continent=="NorthAmerica", continent:="N. America"]
metadata[continent=="SouthAmerica", continent:="S. America"]
#exclude Kenya 2018 outlier
samps<-samps[samps!="SRR14982088"]
#get snp info
a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]


# snpset.cm <- snpgdsLDpruning(genofile,
#                              ld.threshold=0.2,
#                              slide.max.bp = 5000,
#                              autosome.only=FALSE,
#                              sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]],
#                              snp.id=info[chr!="Scaffold_3", snp.id],
#                              maf=3/(2*length(samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]])))
# 
# 
# snpset.id.cm <-unlist(snpset.cm)
# 
# cm.genos<-snpgdsGetGeno(genofile,sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]], snp.id=snpset.id.cm )
# cm.genos[is.na(cm.genos)]<-3
# 
# cm.pca<-prcomp(cm.genos)
# pca.all3<-as.data.table(cm.pca$eigenvect)
# pca.all3[,sample.id:=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]]]
# 
# pca.all3<-merge(pca.all3, metadata[,.(sample.id, Year)])
# 
# set.seed(100)
# dapcTemp_year<-dapc(cm.genos, pca.all3$Year, perc.pca=100, n.da=3)
# ascore_species <- optim.a.score(dapcTemp_year, smart = FALSE, n.sim = 10) 
# 
# dapc_year <- dapc(cm.genos, pca.all3$Year, 
#                   n.pca = 45, n.da = 3)
# scatter.dapc(dapc_year, scree.pca = F, scree.da = F, legend = TRUE, col=friendly_pal("ito_seven"))
# loadingplot(dapc_year$var.contr) 

#what if we include florida in DAPC?

snpset.cm <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"|loc.spec=="FL" , sample.id]],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=3/(2*length(samps[samps%in%metadata[loc.spec=="VA-CM"|loc.spec=="FL", sample.id]])))


snpset.id.cm <-unlist(snpset.cm)

cm.genos<-snpgdsGetGeno(genofile,sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"|loc.spec=="FL", sample.id]], snp.id=snpset.id.cm )
cm.genos[is.na(cm.genos)]<-3

cm.pca<-prcomp(cm.genos)
pca.all3<-as.data.table(cm.pca$eigenvect)
pca.all3[,sample.id:=samps[samps%in%metadata[loc.spec=="VA-CM"|loc.spec=="FL", sample.id]]]

pca.all3<-merge(pca.all3, metadata[,.(sample.id, Year, loc.spec)])
pca.all3[,Year:=as.character(Year)]
pca.all3[loc.spec=="FL", Year:="Florida"]

set.seed(100)
dapcTemp_year<-dapc(cm.genos, pca.all3$Year, perc.pca=100, n.da=3)
ascore_species <- optim.a.score(dapcTemp_year, smart = FALSE, n.sim = 10) 
ascore_species
dapc_year <- dapc(cm.genos, pca.all3$Year, 
                  n.pca = ascore_species$best, n.da = 3)
#scatter.dapc(dapc_year, scree.pca = F, scree.da = F, legend = TRUE, col=friendly_pal("ito_seven"))
#loadingplot(dapc_year$var.contr) 

#functions from https://github.com/JeffWeinell/misc.wrappers/blob/main/R/DAPC_adegenet.R
ggscatter.dapc <- function (x, xax = 1, yax = 2, vartype="df", varname=NULL,axis.title.cex=1, grp = x$grp , cpoint=2, col = adegenet::seasun(length(levels(grp))), txt.leg = levels(grp), label = TRUE, pch = 20, solid = 0.9, hideperimeter=FALSE, show.title=TRUE,scree.da = TRUE, scree.pca = FALSE, posi.da = "bottomright", posi.pca = "bottomleft",bg="white", bg.inset = "white", ratio.da = 0.25, ratio.pca = 0.25, inset.da = 0.02, inset.pca = 0.02, inset.solid = 0.5, onedim.filled = TRUE, mstree = FALSE, lwd = 0.25, lty = 1, segcol = "black", legend = FALSE, posi.leg = "topright", cleg = 1, cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 1, xlim = NULL, ylim = NULL, grid = FALSE, addaxes = TRUE, ltyaxes=2, lwdaxes=0.5, origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, label.inds = NULL, new.pred=NULL){
  if(vartype=="df"){
    ind.vals <- x$ind.coord
    if(is.null(varname)){
      varname  <- "Discriminant function"
    }
  } else {
    if(vartype=="pc"){
      ind.vals <- x$tab
      if(is.null(varname)){
        varname  <- "Principle component"
      }
      mstree   <- FALSE
    } else {
      stop("'vartype' must be either 'df' or 'pc'")
    }
  }
  ### Logical indicating if only one dimension retained
  ONEDIM     <- xax == yax | ncol(ind.vals) == 1
  simple.col <- transp(col[1:length(levels(grp))],solid)
  col        <- transp(rep(col, length(levels(grp))), solid)
  pch        <- rep(pch, length(levels(grp)))
  bg.inset   <- transp(bg.inset, inset.solid)
  ### Posterior assignments of individuals to groups
  if (is.null(grp)) {
    grp <- x$grp
  }
  ### Further checking/updating if there is only one discriminant function
  if (is.null(xax) || is.null(yax)) {
    ## a number followed by L indiciates that the class should be an integer
    xax    <- 1L
    yax    <- ifelse(ncol(ind.vals) == 1L, 1L, 2L)
    ONEDIM <- TRUE
  }
  ### Things to do when more than one PC exists.
  if(!ONEDIM){
    coords.df  <- data.frame(x.coords=ind.vals[, xax],y.coords=ind.vals[, yax],Cluster=grp)
    xlim       <- c(-max(abs(coords.df[,"x.coords"])),max(abs(coords.df[,"x.coords"])))
    ylim       <- c(-max(abs(coords.df[,"y.coords"])),max(abs(coords.df[,"y.coords"])))
    ### Creating columns to hold x and y mean of group that each individual belongs too
    unique.clusters <- levels(grp)
    for(z in unique.clusters){
      rows.temp <- which(coords.df$Cluster==z)
      coords.df[rows.temp,"grp.center.x"] <- mean(coords.df[rows.temp,"x.coords"])
      coords.df[rows.temp,"grp.center.y"] <- mean(coords.df[rows.temp,"y.coords"])
    }
    xy3.df           <- do.call(rbind,lapply(1:nrow(coords.df),FUN=function(z){newpoint(p0=coords.df[z,c("grp.center.x","grp.center.y")], p1=coords.df[z,c("x.coords","y.coords")], c=cstar)}))
    coords.df[,"x3"] <- xy3.df[,1]
    coords.df[,"y3"] <- xy3.df[,2]
    ### blank plotting area
    if(cellipse>0){
      ggscatter.tempA      <- ggplot2::ggplot(coords.df, ggplot2::aes(x=x.coords, y=y.coords,color=Cluster,shape=Cluster,fill=Cluster)) + ggplot2::scale_x_continuous(name=paste(varname,xax)) + ggplot2::scale_y_continuous(name=paste(varname,yax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = bg))  + ggplot2::stat_ellipse(color="white") + ggplot2::theme_classic() + ggplot2::geom_blank()
    } else {
      ggscatter.tempA      <- ggplot2::ggplot(coords.df, ggplot2::aes(x=x.coords, y=y.coords,color=Cluster,shape=Cluster,fill=Cluster)) + ggplot2::scale_x_continuous(name=paste(varname,xax)) + ggplot2::scale_y_continuous(name=paste(varname,yax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = bg)) + ggplot2::theme_classic() + ggplot2::geom_blank()
    }
    # Includes box around plotting area
    ggscatter.tempB      <- ggscatter.tempA + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) 
    # Add reference lines (axes) at x=0 and y=0
    if(addaxes){
      ggscatter.tempB  <- ggscatter.tempB + ggplot2::geom_vline(ggplot2::aes(xintercept=0),linetype=ltyaxes,size=lwdaxes,color="lightgray") + ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype=ltyaxes,size=lwdaxes,color="lightgray")
    }
    # Hide axis ticks and labels
    if(hideperimeter){
      ggscatter.tempC  <- ggscatter.tempB + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
    } else {
      ggscatter.tempC  <- ggscatter.tempB + ggplot2::theme(axis.title.x=ggplot2::element_text(size=(axis.title.cex*12)), axis.text.x=ggplot2::element_text(size=(axis.title.cex*10)), axis.title.y=ggplot2::element_text(size=(axis.title.cex*12)), axis.text.y=ggplot2::element_text(size=(axis.title.cex*10)))
    }
    ### Adding the points. The scale to use for colors of points was defined in ggscatter.tempA so no need to redefine colors here.
    ggscatter.temp0      <- ggscatter.tempC + ggplot2::geom_point(size=cpoint,show.legend=TRUE) + ggplot2::scale_color_manual(values=col) + ggplot2::scale_shape_manual(values=pch) #+ ggplot2::scale_fill_manual(values=col,fill=col)
    # Hide or show legend (guide)
    if(!legend){
      ggscatter.temp1  <- ggscatter.temp0 + ggplot2::theme(legend.position = "none")
    } else {
      ggscatter.temp1  <- ggscatter.temp0 + ggplot2::theme(legend.position = c(0.98,0.98), legend.justification=c("right","top"), legend.background = ggplot2::element_rect(fill="white", size=0.25, linetype="solid")) + ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(fill=simple.col,color="black",size=12,shape=22))) # + ggplot2::scale_fill_manual(values=simple.col) 
    }
    # Add ellipses around clusters
    if(cellipse > 0){
      ggscatter.temp2  <- ggscatter.temp1 + ggplot2::stat_ellipse(level=(cellipse*0.43),type="norm",show.legend=FALSE)
    } else {
      ggscatter.temp2  <- ggscatter.temp1
    }
    # Add 'star' lines from each cluster mean to the coordinates of individuals in the cluster.
    if(cstar > 0){
      ggscatter.temp3 <- ggscatter.temp2 + suppressWarnings(ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = x3, y = y3, xend = grp.center.x, yend = grp.center.y, color = Cluster),show.legend=FALSE))
    } else {
      ggscatter.temp3 <- ggscatter.temp2
    }
    if(mstree){
      meanposi <- apply(x$tab, 2, tapply, grp, mean)
      axes     <- c(xax, yax)
      D        <- dist(meanposi)^2
      tre      <- ade4::mstree(D)
      x0       <- unname(x$grp.coord[tre[, 1], axes[1]])
      y0       <- unname(x$grp.coord[tre[, 1], axes[2]])
      x1       <- unname(x$grp.coord[tre[, 2], axes[1]])
      y1       <- unname(x$grp.coord[tre[, 2], axes[2]])
      tree.df  <- data.frame(xA=x0,yA=y0,xB=x1,yB=y1)
      tree.mat <- cbind(x0,y0,x1,y1)
      coords.df[,"tree.x0"] <- x0
      coords.df[,"tree.y0"] <- y0
      coords.df[,"tree.x1"] <- x1
      coords.df[,"tree.y1"] <- y1
      ggscatter.temp4 <- ggscatter.temp3 + suppressWarnings(ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = tree.x0, y = tree.y0, xend = tree.x1, yend = tree.y1), color = segcol, size=lwd,linetype=lty,show.legend=FALSE))
    } else {
      ggscatter.temp4 <- ggscatter.temp3
    }
    #if(!is.null(label)){
    if(label){
      ggscatter.temp5 <- ggscatter.temp4 + ggplot2::geom_label(data=coords.df,ggplot2::aes(x=grp.center.x,y=grp.center.y,label=Cluster),fill="white",size=(clabel*3),show.legend=FALSE)
    } else {
      ggscatter.temp5 <- ggscatter.temp4
    }
    # return(ggscatter.temp5)
  } else {# If only one PC
    scree.da <- FALSE
    if(ncol(ind.vals) == 1) {
      pcLab <- 1
    } else {
      pcLab <- xax
    }
    ### Data frame with coordinates of individuals and the posterior assignment of individuals to groups (clusters)
    coords.df  <- data.frame(coords=ind.vals[, pcLab],Cluster=grp)
    ### Apply the density function to the individual coordinates for individuals in each group
    ldens <- tapply(X=ind.vals[, pcLab], INDEX=grp, FUN=density)
    ### The actual x (coorinates) and y (density) values that are to be plotted
    allx  <- unlist(lapply(ldens, function(e) e$x))
    ally  <- unlist(lapply(ldens, function(e) e$y))
    ## defining locations for x-axis ticks
    xat0 <- seq(from=round(min(allx)),to=round(max(allx)))
    xat  <- xat0[xat0/2 == round(xat0/2)]
    
    #xpoints <- ind.vals[grp == levels(grp)[i], pcLab]
    #ypoints <- rep(0, sum(grp == levels(grp)[i]))
    
    ### ggplot of PC coordinates of groups.
    if(onedim.filled){
      gg.density.temp  <- ggplot2::ggplot(coords.df, ggplot2::aes(x=coords,color=Cluster,fill=Cluster)) + ggplot2::geom_density() + ggplot2::theme_classic() + ggplot2::scale_color_manual(values=col) + ggplot2::scale_fill_manual(values=col) + ggplot2::scale_x_continuous(breaks=xat,name=paste(varname,xax),limits=range(allx))
    } else {
      gg.density.temp  <- ggplot2::ggplot(coords.df, ggplot2::aes(x=coords,color=Cluster,fill=NA)) + ggplot2::geom_density() + ggplot2::theme_classic() + ggplot2::scale_color_manual(values=col) + ggplot2::scale_fill_manual(values=NA) + ggplot2::scale_x_continuous(breaks=xat,name=paste(varname,xax),limits=range(allx))
    }
    gg.density0       <- gg.density.temp + ggplot2::ylab("Density") + ggplot2::geom_text(ggplot2::aes(x=coords,y=rep(0,length(coords)),label=rep("|",length(coords))),show.legend=FALSE) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1), axis.line=ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size=12), axis.text.y = ggplot2::element_text(size=12))
    ### Add a title
    if(show.title){
      gg.density1 <- gg.density0 + ggplot2::labs(title=paste0("K=",K))
    } else {
      gg.density1 <- gg.density0
    }
    ### Remove legend if legend = FALSE
    if(!legend){
      gg.density2 <- gg.density1 + ggplot2::theme(legend.position = "none")
    } else {
      gg.density2 <- gg.density1 + ggplot2::theme(legend.position = c(0.98,0.98), legend.justification=c("right","top"))
    }
    # Hide axis ticks and labels
    if(hideperimeter){
      gg.density3  <- gg.density2 + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())
    } else {
      gg.density3  <- gg.density2
    }
    
  }
  if(ONEDIM){
    return(gg.density3)
  } else {
    return(ggscatter.temp5)
  }
}
newpoint <-function(p0,p1,c){
  v  <- p1-p0
  d  <- sqrt(sum(v^2))
  u  <- v/d
  d2 <- d*c
  p2 <- p0 + (d2*u)
  p2
}

a<-ggscatter.dapc(x=dapc_year)+scale_color_manual(values=friendly_pal("ito_seven")[c(1:4,6)])+theme_cowplot()+theme(legend.position = "none")
## treemix ###

setwd("/scratch/perickso/private/ind_seq/popgen/treemix/") 

source("treemix_ggplot.R")
theme_set(theme_treemix())

edge0<-read_treemix("bootstrap/5pops500snps_constree_bootrep_8")

b<-plot_treemix(edge0, plot.nodes=F)+  scale_x_continuous(expand = expansion(mult = c(0.1, 0.5))) +theme(axis.title=element_text(size=14, vjust=1), axis.text=element_text(size=10))

pdf("/scratch/perickso/private/ind_seq/Figures/Figure4_DAPC_treemix_noinv_cons.pdf", height=4, width=8)

plot_grid(a,b, nrow=1, labels=c("a", "b"), align="h")
dev.off()

