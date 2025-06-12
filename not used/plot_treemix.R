
setwd("/scratch/perickso/private/ind_seq/popgen/treemix/") # of course this needs to be adjusted
prefix="zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.NorthAmerica.plink.pruneddata.treemix.frq.gz"
library(RColorBrewer)
library(R.utils)
source("plotting_functions.R") # here you need to add the path

par(mfrow=c(2,3))
for(edge in 6:10){
  plot_tree(cex=0.8,paste0(prefix,edge))
  title(paste(edge,"edges"))
}


samps<-data.table(pop=unique(fread("treemix_cluster.txt", header=F)$V3))
samps[,color:=c("#7FFF00", "#7FFF00", "#68228B", "#0000FF", "#EE1289", "#EE1289", "#FF7F00", "#FF7F00", "#0000FF", "#7FFF00")]
write.table(samps, file="treemix_groups.txt", quote=F, row.names=F, col.names=F, sep="\t")
par(mfrow=c(2,3))
for(edge in 0:5){
  plot_resid(stem=paste0(prefix,edge),pop_order="treemix_groups.txt")
}

par(mfrow=c(1,1))

edge=2
plot_tree(paste0(prefix,edge))


par(mfrow=c(2,3))
for(edge in 0:5){
  plot_resid(stem=paste0(prefix,edge),pop_order="treemix_groups.txt")
}


prefix="zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.CM_MIA.plink.pruneddata.treemix.frq.gz"
samps<-data.table(pop=unique(fread("treemix_cluster_CM_MIA.txt", header=F)$V3))
write.table(samps, file="treemix_groups_CM_MIA.txt", quote=F, row.names=F, col.names=F, sep="\t")


par(mfrow=c(2,3))
for(edge in 0:5){
  plot_resid(stem=paste0(prefix,edge),pop_order="treemix_groups_CM_MIA.txt")
}

par(mfrow=c(2,3))
for(edge in 0:5){
  plot_tree(cex=0.8,paste0(prefix,edge))
  title(paste(edge,"edges"))
}

par(mfrow=c(1,1))

edge=4
plot_tree(paste0(prefix,edge))



#try other software from here: https://rdrr.io/github/andrewparkermorgan/popcorn/src/R/treemix.R
library(ggplot2)
source("treemix_ggplot.R")
theme_set(theme_treemix())

edge0<-read_treemix("zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.CM_MIA.plink.pruneddata.treemix.frq.gz0")
edge1<-read_treemix("zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.CM_MIA.plink.pruneddata.treemix.frq.gz1")
edge3<-read_treemix("zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.CM_MIA.plink.pruneddata.treemix.frq.gz3")

edge6<-read_treemix("zaprionus.individual.nosingleton.2023.withID.autosomesonly.noN.NorthAmerica.plink.pruneddata.treemix.frq.gz6")
plot_treemix(edge6)
plot_treemix(edge3, plot.nodes=F)
plot_treemix(edge0)
