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


#read in list of pops as a vector

pops<-fread("/scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt", header=F)$V1

final<-foreach(pop=pops, .errorhandling="remove")%do%{
  
  p.data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/smcpp/", pop, "/", pop, "_cv_10folds_16knots_10000gen_autosomes/", pop, "_CV10fold_final_plot.csv"))
  p.data[,pop:=pop]
  return(p.data)
}
final<-rbindlist(final)
final[,loc:=tstrsplit(pop, split="_")[[1]]]
final[,ky:=x/1000]
final[,type:=ifelse(loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"), "native", "introduced")]
final[loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"),type2:="native"]
final[label=="MIA_2019_June", type2:="introduced-Florida"]
final[grepl("early", label),type2:="introduced-Virginia"]
final[grepl("mid", label),type2:="introduced-Virginia"]

final<-final[!grepl("late", label)]
final[is.na(type2), type2:="introduced"]
final[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
final<-final[order(type2)]
final[,label:=as.factor(label)]
final[,label:=factor(label, levels=rev(levels(final$label)))]

# folds<-foreach(pop=pops)%do%{
#   
#   folds2<-foreach(i=c(0:9), .errorhandling="remove")%do%{
#     p.data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/smcpp/", pop, "/", pop, "_cv_10folds_16knots_10000gen_autosomes/fold", i, "/", pop, "_cv_10folds_16knots_10000gen_fold", i, "_plot.csv"))
#     p.data[,fold:=i]
#     p.data[,pop:=pop]
#     return(p.data)
#     
#   }
#   return(rbindlist(folds2))
# }
# 
# folds<-rbindlist(folds)
# folds[,loc:=tstrsplit(pop, split="_")[[1]]]
# folds[,type:=ifelse(loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"), "native", "introduced")]
# folds[loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"),type2:="native"]
# folds[label=="MIA_2019_June", type2:="introduced-Florida"]
# folds[grepl("early", label),type2:="introduced-early"]
# folds[grepl("mid", label),type2:="introduced-early"]
# 
# folds[grepl("late", label),type2:="introduced-late"]
# folds[is.na(type2), type2:="introduced"]
# folds[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-early",  "native" ))]
# folds<-folds[order(type2)]


smcpp.plot<-ggplot()+
 # geom_line(data=folds[type2%in%c("native", "introduced-Florida", "introduced-early")], aes(x=x, y=y, group =interaction(pop,fold), color=as.factor(type2)), size=0.1)+
  geom_line(data=final[type2%in%c("native", "introduced-Florida", "introduced-Virginia")], aes(x=x, y=y, color=type2, group=label), size=1.5)+  
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +annotation_logticks()+
  labs(x="years before present", y="population size", color=NULL)+
  scale_color_manual(values = friendly_pal("contrast_three"))+
  theme(legend.position = c(0.5, 0.1))



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
metadata[Location=="MIA", Location:="FL"]
metadata[,group:=paste(Location, Year, Season, sep="_")]
groups<-unique(metadata$group)[17:26]

setnames(kinship.table, "ID1", "sample.id")
kinship.table<-merge(kinship.table, metadata, by="sample.id")
setnames(kinship.table, "sample.id", "sample.id1")
setnames(kinship.table, "ID2", "sample.id")
kinship.table<-merge(kinship.table, metadata, by="sample.id", suffixes=c(".pop1", ".pop2"))
kinship.table[,time.diff:=abs(date.pop1-date.pop2)]

#https://www.cell.com/ajhg/pdf/S0002-9297(12)00309-6.pdf

kinship.table[kinship>0.03125&k0<.875, kinship.group:="4th degree"]
kinship.table[kinship>.0625&k0<.75, kinship.group:="3rd degree"]
kinship.table[kinship>.125&k0<.5, kinship.group:="2nd degree"]
kinship.table[kinship>.25&k0<.25, kinship.group:="full sibs"]


kinship.table[is.na(kinship.group), kinship.group:="< 4th degree"]
kinship.table[,kinship.group:=factor(kinship.group, levels=c("full sibs", "2nd degree", "3rd degree", "4th degree", "< 4th degree"))]


k.all<-ggplot(kinship.table[group.pop1%in%groups &group.pop1==group.pop2], aes(x=k0, y=kinship, color=group.pop1))+
  geom_point()+
  labs(color="Population")+
  labs(x="IBD0")+
  scale_color_manual(values = c(friendly_pal("ito_seven"), friendly_pal("bright_seven")))

#kinship over time
k.time<-ggplot(kinship.table[Location.pop1=="VA-CM"&Location.pop2=="VA-CM"], aes(x=time.diff, y=kinship, color=kinship.group))+
  geom_point()+
  labs(x="Days between sampling", color="Relatedness")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  theme(legend.position=c(0.65, 0.8))


kinship.plots<-plot_grid(k.all, k.time, nrow=2, labels=c("B", "C"), align="h", axis="l")

jpeg("/scratch/perickso/private/ind_seq/Figures/Figure3_updated_autosomes.jpeg", height=6, width=10, res=300, units="in")
plot_grid(smcpp.plot, kinship.plots, align="v", axis="b", labels=c("A",""))
dev.off()


#supplemental figure: all plots split 
jpeg("/scratch/perickso/private/ind_seq/Figures/kinship_vs_k0_bypop.jpeg", height=6, width=8, res=300, units="in")
ggplot(kinship.table[group.pop1%in%groups &group.pop1==group.pop2], aes(x=k0, y=kinship, color=kinship.group))+
  geom_point()+
  facet_wrap(~group.pop1, scales="free_x")+
  labs(color="Relatedness")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(limits=c(0,1), breaks=c(0,.5,1))+
  scale_y_continuous(limits=c(0,.25), breaks=c(0.1, 0.2))+
  theme(strip.text.x = element_text(size = 8), axis.text.x=element_text(size=10))+
  labs(x="IBD0")+ theme(panel.spacing = unit(2, "lines"))
dev.off()


