library(data.table)
library(ggplot2)
library(foreach)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", header=T, drop=1)

#first use linux to determine the lowest CV for each chromosome and record here
#grep CV log.admix_unrelated*.[0-9].[0-9].out > CVs_by_chr_unrelated.txt
CV<-fread("/scratch/perickso/private/ind_seq/popgen/CVs_by_chr_unrelated.txt")
CV[,chr:=tstrsplit(V1, split="[.]")[[3]]]
CV[,K:=tstrsplit(V1, split="[.]")[[4]]]

CV.min<-CV[,.(min.CV.k=K[V4==min(V4)]), .(chr)]



#reorder k groups so that africa is always k1

y<-foreach(chrom=c(1:5))%do%{
  min.k<-CV.min[chr==chrom,min.CV.k]
  if(chrom==3){
    names=fread("/scratch/perickso/private/ind_seq/popgen/admix_unrelated_bychr.3.plink.fam", header=F)$V2
  }else{
  names=fread("/scratch/perickso/private/ind_seq/popgen/admix_unrelated_bychr.1.plink.fam", header=F)$V2
  }
  admix<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/admix_unrelated_bychr.",chrom, ".pruneddata.", min.k, ".Q" ), header=F)
  names(admix)=paste0("k", c(1:min.k))
  admix[,sample.id:=names]
  admix.melt<-melt(admix, id.vars="sample.id", value.name="proportion",variable.name="k.group")
  admix.melt[,k.total:=min.k]
  admix.melt<-merge(admix.melt, metadata, by="sample.id")
  k.props<-admix.melt[continent=='Africa', .(avg.prop=mean(proportion)), .(k.group)]
  top.k<-k.props[order(avg.prop, decreasing=T)]
  top.k[,new.k.group:=1:nrow(top.k)]
  admix.melt<-merge(admix.melt, top.k, by="k.group")
  admix.melt[,chr:=chrom]
  return(admix.melt)
}
y<-rbindlist(y)


y[continent=="Africa", loc.spec:="Africa"]
y[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
y[,k.label:=paste0("chr ", chr, ", k=", k.total)]

fwrite(y, file="/scratch/perickso/private/ind_seq/popgen/admixture_data_for_figure.csv")

#for paper figure:

jpeg('/scratch/perickso/private/ind_seq/Figures/admixture_all_by_chrom_unrelated.jpg', height=6, width=10, units="in", res=2400)

ggplot(y[loc.spec%in%c("Africa", "FL", "VA-HPO", "VA-CM")], aes(fill=as.factor(new.k.group), y=proportion, x=sample.id)) + 
  geom_bar(position="fill", stat="identity")+
  facet_grid(k.label~loc.spec, scales = "free_x", space = "free_x")+ 
  theme(axis.text.x = element_blank())+
  guides(fill="none")+
  labs(x=NULL)+
  # theme(strip.text.x = element_text(size=6,angle=75))+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())

dev.off()

#for 