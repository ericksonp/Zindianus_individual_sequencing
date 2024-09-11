#plot scaffold 3 data from different programs as well as whole genome comparison

library(data.table)
library(ggplot2)
library(ggpubfigs)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(ggbio)
#library(gggenomes)
#everything will focus on 3 "populations": AFrica, Florida, and Virginia (combined females)
#PIXY: Dxy and Pi

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
pixy.fst

#fwrite(pixy.fst[comparison=="FL-VA"], file='/scratch/perickso_shared/alexandra/FST_VAvsFL_windowed.csv') #save for Alexandra
#add in high IHS points at the top

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[order(IHS)] #973443

ihs.high<-ihs.VA[IHS>5&CHR=="Scaffold_3"&POSITION<1500000, POSITION]
snps.to.plot<-data.table(x=ihs.high,
                         y=0.8)

a<-ggplot()+
  geom_line(data=pixy.fst[pop1!="other"&pop2!="other"][chromosome=="Scaffold_3"][window.mid>250000&window.mid<1100000], aes(x=window.mid, 
                y=avg_wc_fst, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("F"[ST]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  geom_point(data=snps.to.plot, aes(x=x, y=y), alpha=0.25, stroke=NA)


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


b<-ggplot(pixy.dxy[pop1!="other"&pop2!="other"][chromosome=="Scaffold_3"][window.mid>250000&window.mid<1100000])+
  geom_line(aes(x=window.mid,
                y=avg_dxy, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("D"[xy]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  guides(color="none")

  





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
pixy.pi


c<-ggplot(pixy.pi[pop!="other"][chromosome=="Scaffold_3"][window.mid>250000&window.mid<1100000])+
  geom_line(aes(x=window.mid, 
                y=avg_pi, 
                color=pop, 
                group=pop))+
  labs(x=NULL, 
       y="π",
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  theme(axis.text.x=element_blank())

#tajima's D
thetas<-foreach(pop=c("Africa", "FL", "VA"), .errorhandling="remove")%do%{
  theta<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/angsd/", pop, ".thetasWindow5000_5000.pestPG"))
  theta[,group:=pop]
  return(theta)
}

thetas<-rbindlist(thetas)

thetas[,chrom:=as.numeric(tstrsplit(Chr, split="_")[[2]])]




d<-ggplot(thetas[chrom==3][WinCenter>250000&WinCenter<1100000])+
  geom_line(aes(x=WinCenter, y=Tajima, color=group, group=group))+
  labs(x=NULL, 
       y="Tajima's D", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  geom_hline(yintercept = 0, color="grey50", linetype="dashed")+
  theme(axis.text.x=element_blank())+
  guides(color="none")


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


e<-ggplot(depth.sum[chr==3&win>250000&win<1100000&label%in%(c("Africa", "FL", "VA"))], 
       (aes(x=win, 
            y=mean.depth,
            ymin=mean.depth-sd.depth,
            ymax=mean.depth+sd.depth,
            group=label, 
            color=label)))+
  geom_line() +
  labs(x="Position", 
       y="mean relative \nsequencing depth", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE))+
  guides(color="none")




jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_scaffold_3.jpg", 
     height=10, 
     width=8,
     units="in", 
     res=300 )
plot_grid(a,b,c,d,e, nrow=5, labels=c("A", "B", "C", "D", "E"), align = 'v', axis = 'lr')
dev.off()

#exploratory: try to get gene annotations! 
# gff<-"/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff"
# 
# ann<-read_feats(gff)
# ann.dt<-as.data.table(ann)
# 
# f<-gggenomes(ann.dt[seq_id=="Scaffold_3"&start>250000&start<1000000])+geom_gene(size=10)+ geom_seq_label()+lims(x=c(0,750000))
# 
# 
# plot_grid(a,b,c,d,e,f, nrow=6, labels=c("A", "B", "C", "D", "E", "F"), align = 'v', axis = 'lr')
