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

pixy.fst.sum<-pixy.fst[chromosome!="Scaffold_3"&pop1!="other"&pop2!="other",.(q95=quantile(avg_wc_fst, 0.95, na.rm=T), q99=quantile(avg_wc_fst, 0.99, na.rm=T)), .(comparison)]

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-wgscan.ihs$ihs
ihs.high2<-ihs.VA[CHR=="Scaffold_2"&POSITION>26590000&POSITION<26700000, POSITION]
ihs.high5<-ihs.VA[CHR=="Scaffold_5"&POSITION>7840000&POSITION<7920000, POSITION]


ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons2<-ann[chr=="Scaffold_2"&info=="exon"&((start>min(ihs.high2)&start<max(ihs.high2))|(end>min(ihs.high2)&end<max(ihs.high2)))]
exons5<-ann[chr=="Scaffold_5"&info=="exon"&((start>min(ihs.high5)&start<max(ihs.high5))|(end>min(ihs.high5)&end<max(ihs.high5)))]
genes2<-ann[chr=="Scaffold_2"&info=="gene"&((start>min(ihs.high2)&start<max(ihs.high2))|(end>min(ihs.high2)&end<max(ihs.high2)))]
genes5<-ann[chr=="Scaffold_5"&info=="gene"&((start>min(ihs.high5)&start<max(ihs.high5))|(end>min(ihs.high5)&end<max(ihs.high5)))]



a<-ggplot()+
  geom_line(data=pixy.fst[pop1!="other"&pop2!="other"][chromosome=="Scaffold_2"][window.mid>26590000&window.mid<26700000], aes(x=window.mid, 
                y=avg_wc_fst, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("F"[ST]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  geom_segment(data=genes2, aes(x=start, xend=end, y=.8, yend=.8))+
  geom_rect(data=exons2, aes(xmin=start,xmax=end, ymin=.7, ymax=.9, fill=ifelse(gene_name=="ANN06929", "black", "grey")))+
  scale_fill_manual(values = c(friendly_pal("ito_seven")[2], "grey40"))+
  guides(color="none", fill="none")+
  geom_hline(data=pixy.fst.sum, aes(yintercept=q99, color=comparison), linetype="dashed")
  
  
b<-ggplot()+
  geom_line(data=pixy.fst[pop1!="other"&pop2!="other"][chromosome=="Scaffold_5"][window.mid>7800000&window.mid<7950000], 
            aes(x=window.mid, 
                y=avg_wc_fst, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("F"[ST]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  geom_segment(data=genes5, aes(x=start, xend=end, y=.8, yend=.8))+
  geom_rect(data=exons5[gene_name!="ANN03542"], aes(xmin=start,xmax=end, ymin=.7, ymax=.9, fill=gene_name))+
  #plot other gene in gray
  geom_rect(data=exons5[gene_name=="ANN03542"], aes(xmin=start,xmax=end, ymin=.7, ymax=.9), fill="grey50")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  guides(fill="none")+  
  geom_hline(data=pixy.fst.sum, aes(yintercept=q99, color=comparison), linetype="dashed")


  

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
pixy.dxy.sum<-pixy.dxy[chromosome!="Scaffold_3"&pop1!="other"&pop2!="other",.(q95=quantile(avg_dxy, 0.95, na.rm=T), q99=quantile(avg_dxy, 0.99, na.rm=T), median=median(avg_dxy, na.rm=T)), .(comparison)]


c<-ggplot(pixy.dxy[pop1!="other"&pop2!="other"][chromosome=="Scaffold_2"][window.mid>26590000&window.mid<26700000])+
  geom_line(aes(x=window.mid,
                y=avg_dxy, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("D"[xy]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,.025))+
  geom_hline(data=pixy.dxy.sum, aes(yintercept=q99, color=comparison), linetype="dashed")


d<-ggplot(pixy.dxy[pop1!="other"&pop2!="other"][chromosome=="Scaffold_5"][window.mid>7800000&window.mid<7950000])+
  geom_line(aes(x=window.mid,
                y=avg_dxy, 
                color=comparison, 
                group=comparison))+
  labs(x=NULL, 
       y=expression("D"[xy]))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,0.025))+
  geom_hline(data=pixy.dxy.sum, aes(yintercept=q99, color=comparison), linetype="dashed")




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

pixy.pi.sum<-pixy.pi[chromosome!="Scaffold_3"&pop!="other",.(q05=quantile(avg_pi, 0.05, na.rm=T), q01=quantile(avg_pi, 0.1, na.rm=T)), .(pop)]


e<-ggplot(pixy.pi[pop!="other"][chromosome=="Scaffold_2"][window.mid>26590000&window.mid<26700000])+
  geom_line(aes(x=window.mid, 
                y=avg_pi, 
                color=pop, 
                group=pop))+
  labs(x=NULL, 
       y="π",
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,0.025))+
  geom_hline(data=pixy.pi.sum, aes(yintercept=q01, color=pop), linetype="dashed")



f<-ggplot(pixy.pi[pop!="other"][chromosome=="Scaffold_5"][window.mid>7800000&window.mid<7950000])+
  geom_line(aes(x=window.mid, 
                y=avg_pi, 
                color=pop, 
                group=pop))+
  labs(x=NULL, 
       y="π",
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  theme(axis.text.x=element_blank())+
  lims(y=c(0,0.025))+
  geom_hline(data=pixy.pi.sum, aes(yintercept=q01, color=pop), linetype="dashed")


#tajima's D
thetas<-foreach(pop=c("Africa", "FL", "VA"), .errorhandling="remove")%do%{
  theta<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/angsd/", pop, ".thetasWindow5000_5000.pestPG"))
  theta[,group:=pop]
  return(theta)
}

thetas<-rbindlist(thetas)

thetas[,chrom:=as.numeric(tstrsplit(Chr, split="_")[[2]])]

thetas.sum<-thetas[chrom!="Scaffold_3",.(q05=quantile(Tajima, 0.05, na.rm=T), q01=quantile(Tajima, 0.1, na.rm=T)), .(group)]



g<-ggplot(thetas[chrom==2][WinCenter>26590000&WinCenter<26700000])+
  geom_line(aes(x=WinCenter, y=Tajima, color=group, group=group))+
  labs(x=NULL, 
       y="Tajima's D", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  #geom_hline(yintercept = 0, color="grey50", linetype="dashed")+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(-2.5, 2))+
  geom_hline(data=thetas.sum, aes(yintercept=q01, color=group), linetype="dashed")


h<-ggplot(thetas[chrom==5][WinCenter>7800000&WinCenter<7950000])+
  geom_line(aes(x=WinCenter, y=Tajima, color=group, group=group))+
  labs(x=NULL, 
       y="Tajima's D", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6])+
  #geom_hline(yintercept = 0, color="grey50", linetype="dashed")+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(-2.5, 2))+
  geom_hline(data=thetas.sum, aes(yintercept=q01, color=group), linetype="dashed")



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


i<-ggplot(depth.sum[chr==2&win>26590000&win<26700000&label%in%(c("Africa", "FL", "VA"))], 
       (aes(x=win/1000000, 
            y=mean.depth,
            #ymin=mean.depth-sd.depth,
            #ymax=mean.depth+sd.depth,
            group=label, 
            color=label)))+
  geom_line() +
  labs(x="Chr. 2 Position (Mb)", 
       y="mean relative \nsequencing depth", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE))+
  guides(color="none")+
  lims(y=c(0.5,2.0))



j<-ggplot(depth.sum[chr==5&win>7800000&win<7950000&label%in%(c("Africa", "FL", "VA"))], 
          (aes(x=win/1000000, 
               y=mean.depth,
               #ymin=mean.depth-sd.depth,
               #ymax=mean.depth+sd.depth,
               group=label, 
               color=label)))+
  geom_line() +
  labs(x="Chr. 5 Position (Mb)", 
       y="mean relative \nsequencing depth", 
       color="Population")+
  scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE))+
  guides(color="none")+
  lims(y=c(0.5,2.0))

left<-plot_grid(a,c,e,g,i, nrow=5, labels=c("a", "c", "e", "g", "i"), align="v", axis = "lr")
right<-plot_grid(b,d,f,h,j, nrow=5, labels=c("b", "d", "f", "h", "j"), align="v", axis = "lr")

#jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_scaffold_2and5.jpg", height=8, width=8, units="in", res=300 )

pdf("/scratch/perickso/private/ind_seq/popgen/plots/Figure7_diversity_stats_scaffold_2and5.pdf", height=8, width=10)
plot_grid(left, right, nrow=1, rel_widths=c(0.42, 0.58), align="h", axis="tb")
dev.off()
