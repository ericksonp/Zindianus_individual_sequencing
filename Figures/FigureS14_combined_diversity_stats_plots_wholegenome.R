#plot scaffold 3 data from different programs as well as whole genome comparison

library(data.table)
library(ggplot2)
library(ggpubfigs)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)

#everything will focus on 3 "populations": AFrica, Florida, and Virginia (combined females)
#PIXY: Dxy and Pi

pixy.dxy<-foreach(i=c(1:5))%do%{
  if(i==3) {
    data<-fread("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_3_females_dxy.txt")
  } else {
    data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_", i, "_dxy.txt"))
  }
  return(data)
}
pixy.dxy<-rbindlist(pixy.dxy)
pixy.dxy[,window.mid:=(window_pos_1+window_pos_2)/2]
pixy.dxy[,index:=rleid(chromosome, window_pos_1)]
pixy.dxy[,comparison:=paste(pop1, pop2, sep="-")]
pixy.dxy.rect<-data.table(starts=c(min(pixy.dxy[chromosome=="Scaffold_2", index]),
                                   min(pixy.dxy[chromosome=="Scaffold_4", index])),
                          ends=c(max(pixy.dxy[chromosome=="Scaffold_2", index]),
                                 max(pixy.dxy[chromosome=="Scaffold_4", index])))



a<-ggplot()+
  geom_rect(data=pixy.dxy.rect, aes(xmin=starts, xmax=ends, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_line(data=pixy.dxy[pop1!="other"&pop2!="other"],aes(x=index,
                y=avg_dxy, 
                color=comparison, 
                group=comparison), linewidth=0.2)+
  lims(y=c(0,.06))+
  labs(x=NULL, 
       y=expression("D"[xy]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.x=element_blank())
  
pixy.fst<-foreach(i=c(1:5))%do%{
  if(i==3) {
    data<-fread("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_3_females_fst.txt")
  } else {
    data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_", i, "_fst.txt"))
  }
  return(data)
}
pixy.fst<-rbindlist(pixy.fst)
pixy.fst[,window.mid:=(window_pos_1+window_pos_2)/2]
pixy.fst[,index:=rleid(chromosome, window_pos_1)]
pixy.fst[,comparison:=paste(pop1, pop2, sep="-")]
pixy.fst.rect<-data.table(starts=c(min(pixy.fst[chromosome=="Scaffold_2", index]),
                                   min(pixy.fst[chromosome=="Scaffold_4", index])),
                          ends=c(max(pixy.fst[chromosome=="Scaffold_2", index]),
                                 max(pixy.fst[chromosome=="Scaffold_4", index])))




b<-ggplot()+
  geom_rect(data=pixy.fst.rect, aes(xmin=starts, xmax=ends, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_line(data=pixy.fst[pop1!="other"&pop2!="other"], 
            aes(x=index, 
                y=avg_wc_fst, 
                color=comparison, 
                group=comparison), linewidth=0.2)+
  labs(x=NULL, 
     y=expression("F"[ST]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.x=element_blank())



pixy.pi<-foreach(i=c(1:5))%do%{
  if(i==3) {
    data<-fread("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_3_females_pi.txt")
  } else {
    data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_", i, "_pi.txt"))
  }
  return(data)
}
pixy.pi<-rbindlist(pixy.pi)
pixy.pi[,window.mid:=(window_pos_1+window_pos_2)/2]
pixy.pi[,index:=rleid(chromosome, window_pos_1)]
pixy.pi.rect<-data.table(starts=c(min(pixy.pi[chromosome=="Scaffold_2", index]),
                                   min(pixy.pi[chromosome=="Scaffold_4", index])),
                          ends=c(max(pixy.pi[chromosome=="Scaffold_2", index]),
                                 max(pixy.pi[chromosome=="Scaffold_4", index])))


c<-ggplot()+
  geom_rect(data=pixy.pi.rect, aes(xmin=starts, xmax=ends, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_line(data=pixy.pi[pop!="other"],
            aes(x=index, 
                y=avg_pi, 
                color=pop, 
                group=pop), linewidth=0.2)+
  lims(y=c(0,.06))+
  labs(x=NULL, 
       y="π",
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  theme(axis.text.x=element_blank())+
  theme(axis.text.x=element_blank())


#tajima's D
thetas<-foreach(pop=c("Africa", "FL", "VA"), .errorhandling="remove")%do%{
  theta<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/angsd/", pop, ".thetasWindow5000_5000.pestPG"))
  theta[,group:=pop]
  return(theta)
}
thetas<-rbindlist(thetas)
thetas[,chrom:=as.numeric(tstrsplit(Chr, split="_")[[2]])]
thetas<-thetas[chrom<=5]
thetas<-thetas[order(WinCenter)][order(Chr)]
thetas[,index:=rleid(Chr, WinCenter)]
thetas.rect<-data.table(starts=c(min(thetas[Chr=="Scaffold_2", index]),
                                  min(thetas[Chr=="Scaffold_4", index])),
                         ends=c(max(thetas[Chr=="Scaffold_2", index]),
                                max(thetas[Chr=="Scaffold_4", index])))


d<-ggplot()+
  geom_rect(data=thetas.rect, aes(xmin=starts, xmax=ends, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_line(data=thetas, aes(x=index, 
                             y=Tajima, 
                             color=group, 
                             group=group), linewidth=0.2)+
  labs(x=NULL, 
       y="Tajima's D", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  geom_hline(yintercept = 0, color="grey50", linetype="dashed")+
  theme(axis.text.x=element_blank())+
  theme(axis.text.x=element_blank())



#depth
depth<-fread("/scratch/perickso/private/ind_seq/sv/depth_analysis_5kb.csv")
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", header=T, drop=1)

depth<-merge(depth, metadata, by="sample.id")

#look at summarized data

depth[,label:=loc.spec]
depth[label=="VA-HPO"|label=="VA-CM", label:="VA"]
depth.sum<-depth[,.(mean.depth=mean(rel_depth),
                    sd.depth=sd(rel_depth)),
                 .(label, chr, win)]

depth.sum<-depth.sum[order(win)][order(chr)]
depth.sum[,index:=rleid(win, chr)]

depth.rect<-data.table(starts=c(min(depth.sum[chr==2, index]),
                                min(depth.sum[chr==4, index])),
                       ends=c(max(depth.sum[chr==2, index]),
                              max(depth.sum[chr==4, index])))

e<-ggplot()+ 
  geom_rect(data=depth.rect, aes(xmin=starts, xmax=ends, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_line(data=depth.sum[label%in%(c("Africa", "FL", "VA"))], 
            aes(x=index, 
                 y=mean.depth,
                 group=label, 
                 color=label), linewidth=0.2) +
  lims(y=c(0,2))+
  labs(x="Genome position", 
       y="mean relative \nsequencing depth", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  theme(axis.text.x=element_blank())

jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_whole_genome.jpg", 
     height=10, 
     width=8,
     units="in", 
     res=300 )
plot_grid(a,b,c,d,e, nrow=5, labels=c("A", "B", "C", "D", "E"), align = 'v', axis = 'lr')
dev.off()


