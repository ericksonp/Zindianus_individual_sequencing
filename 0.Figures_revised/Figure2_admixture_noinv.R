library(data.table)
library(ggplot2)
library(foreach)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)

#first look at whole genome data, not removing related individuals
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", header=T, drop=1)

#first use linux to determine the lowest CV for each chromosome and record here
#grep CV log.admix_autosomes_noinv..[0-9].out > CVs_admix_autosomes_noinv.txt
CV<-fread("/scratch/perickso/private/ind_seq/popgen/CVs_admix_autosomes_noinv.txt")
CV[,chr:=tstrsplit(V1, split="[.]")[[3]]]
CV[,K:=tstrsplit(V1, split="[.]")[[4]]]

CV.min<-CV[,.(min.CV.k=K[V4==min(V4)]), .(chr)]



#reorder k groups so that africa is always k1
names=fread("/scratch/perickso/private/ind_seq/popgen/admix_autosomes_noinv.plink.fam", header=F)$V2
y<-foreach(k=c(2:6))%do%{

  admix<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/admix_autosomes_noinv.pruneddata.", k, ".Q" ), header=F)
  names(admix)=paste0("k", c(1:k))
  admix[,sample.id:=names]
  admix.melt<-melt(admix, id.vars="sample.id", value.name="proportion",variable.name="k.group")
  admix.melt[,k.total:=k]
  admix.melt<-merge(admix.melt, metadata, by="sample.id")
  k.props<-admix.melt[continent=='Africa', .(avg.prop=mean(proportion)), .(k.group)]
  top.k<-k.props[order(avg.prop, decreasing=T)]
  top.k[,new.k.group:=1:nrow(top.k)]
  admix.melt<-merge(admix.melt, top.k, by="k.group")
  return(admix.melt)
}
y<-rbindlist(y)


y[continent=="Africa", loc.spec:="Africa"]
y[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
y[,k.label:=paste0( "k = ", k.total)]

fwrite(y, file="/scratch/perickso/private/ind_seq/popgen/admixture_data_for_figure_autosomes_noinv_mac3.csv")

#for paper figure:
y[,loc.label:=loc.spec]
y[loc.label=="VA-CM"|loc.spec=="VA-HPO", loc.label:=paste0(loc.spec, "-", Year )]
y[,loc.label:=factor(loc.label, levels=c("Africa", "FL", "VA-CM-2017", "VA-CM-2018", "VA-CM-2019", "VA-CM-2020", "VA-HPO-2019", "VA-HPO-2020"))]

pdf('/scratch/perickso/private/ind_seq/Figures/Figure2_admixture_autosomes_noinv.pdf', height=5, width=8)

#jpeg('/scratch/perickso/private/ind_seq/Figures/admixture_autosomes_noinv.jpg', height=6, width=10, units="in", res=2400)

ggplot(y[loc.spec%in%c("Africa", "FL", "VA-HPO", "VA-CM")], aes(fill=as.factor(new.k.group), y=proportion, x=sample.id)) + 
  geom_bar(position="fill", stat="identity")+
  facet_grid(k.label~loc.label, scales = "free_x", space = "free_x")+ 
  theme(axis.text.x = element_blank())+
  guides(fill="none")+
  labs(x=NULL)+
  # theme(strip.text.x = element_text(size=6,angle=75))+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(),
        strip.text.x = element_text(size = 8))
dev.off()
