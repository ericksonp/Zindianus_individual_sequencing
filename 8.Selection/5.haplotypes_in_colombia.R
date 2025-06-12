# are haplotypes in colombia?
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)
library(gdsfmt)
library(SNPRelate)
library(lattice)
library(tidyr)
library(stringr)
library(lubridate)
library(viridis)
library(ggpubfigs)
library(rehh)
library(scales)

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[CHR=="Scaffold_5"][order(IHS, decreasing=T)] # 7875359

ihs.high5<-ihs.VA[CHR=="Scaffold_5"&POSITION>7840000&POSITION<7920000, POSITION]

load("/scratch/perickso/private/ind_seq/popgen/scaffold_5haplotype_table.Rdat")

haplos.to.plot5<-haplos.melt[pos%in%ihs.high5]

focal.geno5<-haplos.to.plot5[pos==7875359] #ordering based on IHS peak
focal.geno5[,fixed.geno:=genotype]
haplos.to.plot5<-merge(haplos.to.plot5, focal.geno5[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot5<-haplos.to.plot5[loc.spec%in%c("Africa","FL", "Colombia", "HI")][order(loc.spec, -fixed.geno)]


haplos.to.plot5[,haplo.index:=rleid(haplo.id)]
haplos.to.plot5[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot5<-unique(haplos.to.plot5$pos)


haplo.plot5<-ggplot(haplos.to.plot5)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+
  
  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(strip.text.y = element_text(size = 8))

#colombia and Hawaii mostly have the derived haplotype for scaffold 5



ihs.high2<-ihs.VA[CHR=="Scaffold_2"&POSITION>26590000&POSITION<26700000, POSITION]
load("/scratch/perickso/private/ind_seq/popgen/scaffold_2haplotype_table.Rdat")


#haplos.to.plot<-haplos.melt[pos%in%fst.high|pos%in%ihs.high]
haplos.to.plot2<-haplos.melt[pos%in%ihs.high2]

focal.geno2<-haplos.to.plot2[pos==26609601] #ordering based on IHS peak
focal.geno2[,fixed.geno:=genotype]
haplos.to.plot2<-merge(haplos.to.plot2, focal.geno2[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot2<-haplos.to.plot2[loc.spec%in%c("Africa", "HI", "Colombia", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot2[,haplo.index:=rleid(haplo.id)]
haplos.to.plot2[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot2<-unique(haplos.to.plot2$pos)


haplo.plot2<-ggplot(haplos.to.plot2)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+
  
  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(strip.text.y = element_text(size = 8))

haplo.plot2

#virginia haplotype is present in colombia and hawaii but not fixed
