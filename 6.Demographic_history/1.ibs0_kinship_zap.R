library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(lattice)
library(tidyr)
library(stringr)
library(lubridate)
library(viridis)
library(ggpubfigs)

#IBSO/kinship
geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gds" , allow.fork=T)


a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
autosomal.snps<-info[chr!="Scaffold_3", snp.id]
ibd<-snpgdsIBDMoM(geno, autosome.only=F, missing.rate=0.10, maf=0.01, snp.id=autosomal.snps)
#ibd<-snpgdsIBDKING(geno,type="KING-homo", autosome.only=F, missing.rate=0.01)
#ibd.king.robust<-snpgdsIBDKING(geno,type="KING-robust", autosome.only=F)

kinship.table<-as.data.table(snpgdsIBDSelection(ibd))

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
metadata[,date:=mdy(date)]
metadata[,group:=paste(Location, Year, Season, sep="_")]
groups<-unique(metadata$group)[17:26]

setnames(kinship.table, "ID1", "sample.id")
kinship.table<-merge(kinship.table, metadata, by="sample.id")
setnames(kinship.table, "sample.id", "sample.id1")
setnames(kinship.table, "ID2", "sample.id")
kinship.table<-merge(kinship.table, metadata, by="sample.id", suffixes=c(".pop1", ".pop2"))
kinship.table[,time.diff:=abs(date.pop1-date.pop2)]

kinship.table[kinship>.0625&k0<.75&group.pop1!=group.pop2, kinship.group:="3rd-degree"]
kinship.table[kinship>.0625&k0<.75&group.pop1==group.pop2, kinship.group:="3rd-degree"]
kinship.table[kinship>.125&k0<.5&group.pop1!=group.pop2, kinship.group:="2nd-degree"]
kinship.table[kinship>.125&k0<.5&group.pop1==group.pop2, kinship.group:="2nd-degree"]
kinship.table[kinship>.25&k0<.25&group.pop1==group.pop2, kinship.group:="full sibs"]

kinship.table[is.na(kinship.group), kinship.group:="not related"]
#plot all pairwise-comparisons

#just within-population comparisons
jpeg("/scratch/perickso/private/ind_seq/Figures/kinship_vs_k0_bypop.jpeg", height=6, width=8, res=300, units="in")
ggplot(kinship.table[group.pop1%in%groups &group.pop1==group.pop2], aes(x=k0, y=kinship, color=kinship.group))+
    geom_point(alpha=0.25)+
    facet_wrap(~group.pop1)+
  labs(color="Relatedness")+
  scale_color_manual(values = friendly_pal("bright_seven"))
dev.off()

#kinship over time
jpeg("/scratch/perickso/private/ind_seq/Figures/kinship_over_time.jpeg", height=5, width=5, res=300, units="in")
ggplot(kinship.table[Location.pop1=="VA-CM"&Location.pop2=="VA-CM"], aes(x=time.diff, y=kinship, color=kinship.group))+
  geom_point()+
  labs(x="Days between sampling", color="Relatedness")+
  scale_color_manual(values = friendly_pal("bright_seven"))

dev.off()



jpeg("/scratch/perickso/private/ind_seq/Figures/kinship_over_time_split_by_year.jpeg", height=5, width=5, res=300, units="in")
ggplot(kinship.table[Location.pop1=="VA-CM"&Location.pop2=="VA-CM"&Year.pop1==Year.pop2], aes(x=time.diff, y=kinship, color=kinship.group))+
  geom_point()+
  labs(x="Days between sampling", color="Relatedness")+ 
  scale_color_manual(values = friendly_pal("bright_seven"))
dev.off()


#combine last two into a single plot


a<-ggplot(kinship.table[Location.pop1=="VA-CM"&Location.pop2=="VA-CM"], aes(x=time.diff, y=kinship, color=kinship.group))+
  geom_point()+
  labs(x="Days between sampling", color="Relatedness")+
  scale_color_manual(values = friendly_pal("bright_seven"))+
  guides(color="none")
b<-ggplot(kinship.table[Location.pop1=="VA-CM"&Location.pop2=="VA-CM"&Year.pop1==Year.pop2], aes(x=time.diff, y=kinship, color=kinship.group))+
  geom_point()+
  labs(x="Days between sampling", color="Relatedness")+ 
  scale_color_manual(values = friendly_pal("bright_seven"))+
  facet_wrap(~Year.pop1)

plot_grid(a,b, nrow=1, labels=c("A", "B"), rel_widths=c(0.4, 0.6))
