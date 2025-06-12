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

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#get file with windows containing high ihs, fst, and bayepass scores
snptags<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_snpstoplot.csv")


#FST
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")

fst.plot<-ggplot()+
  geom_rect(data=snptags, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=z, aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y=expression("Florida-Virginia F"[ST]))+
  theme(axis.text.x=element_blank())+  guides(color="none")+
  scale_y_continuous(limits=c(0,1), breaks=c(0, .25, .5, .75))



#Bayescan

bp.va.fl<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/updated_data_for_manhattan.txt")
bp.va.fl[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(old.snp.id=a$snp.id,
                 CHR=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,Scaffold:=as.integer(tstrsplit(CHR, split="_")[[2]])]
bp.va.fl<-merge(bp.va.fl, info, by=c("pos", "Scaffold"))

bp.plot<-ggplot()+  
  geom_rect(data=snptags, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=bp.va.fl, aes(x=old.snp.id, y=M_XtX, color=as.factor(Scaffold)))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y="XtX")+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,9.5))

#IHS

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VAall<-as.data.table(wgscan.ihs$ihs)

setnames(info, "pos", "POSITION")
ihs.VAall<-merge(ihs.VAall, info, by=c("CHR", "POSITION"))


ihs.plot<-ggplot()+
  geom_rect(data=snptags, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=ihs.VAall,aes(x=old.snp.id, y=IHS, color=CHR))+
  scale_color_manual(values = friendly_pal("ito_seven"), labels=c("1", "2", "3", "4", "5"))+
  labs(x="SNP #", y="IHS", color="Chr.")+
  scale_x_continuous(label=scientific)+
  theme(legend.position = "bottom")+
  lims(y=c(-4,8))
  

left<-plot_grid(fst.plot, bp.plot, ihs.plot, nrow=3, align="v", axis="lr", rel_heights =c(0.3,0.3,0.4), labels=c("A", "B", "C"))


#middle of plot-zoomed in images

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_FL.Rdat")
ihs.FL<-as.data.table(wgscan.ihs$ihs)
ihs.FL[,pop:="Florida"]

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_NorthAmerica.Rdat")
ihs.NA<-as.data.table(wgscan.ihs$ihs)
ihs.NA[,pop:="North America"]

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_Africa.Rdat")
ihs.AF<-as.data.table(wgscan.ihs$ihs)
ihs.AF[,pop:="Africa"]

ihs<-rbindlist(list(ihs.FL, ihs.NA, ihs.VA, ihs.AF))
ihs[, scaffold:=tstrsplit(CHR, split="__")[[1]]]
ihs[,abs.IHS:=abs(IHS)]
ihs<-ihs[order(POSITION)][order(scaffold)]
ihs[,snp.id:=c(1:nrow(ihs))]
ihs[,pop:=factor(pop, levels=c("North America", "Virginia", "Florida", "Africa"))]
ihs<-ihs[order(pop)]


fst.sc3<-ggplot(z[chromosome=="Scaffold_3"& position<2000000])+geom_point(aes(x=position, y=fst.snp))+
  #scale_color_manual(values = friendly_pal("ito_seven")[3])+
  labs(x=NULL, y=NULL)+
         theme(axis.text.x=element_blank())+
  guides(color="none")+
  scale_y_continuous(limits=c(0,1), breaks=c(0, .25, .5, .75))
  

bp.sc3<-ggplot(bp.va.fl[CHR=="Scaffold_3"& pos<2000000])+
  geom_point(aes(x=pos, y=M_XtX))+
  #scale_color_manual(values = friendly_pal("ito_seven")[3])+
  labs(x=NULL, y=NULL)+
  #theme(axis.text.x=element_blank())+
  guides(color="none")+
  theme(axis.text.x=element_blank())+
  lims(y=c(0,9.5))
  

  


ihs.sc3<-ggplot()+
  geom_point(data=ihs[scaffold=="Scaffold_3"&POSITION<2000000&pop!="North America"], aes(x=POSITION, y=IHS, color=pop), alpha=0.5)+
  #geom_point(color="white", pch=21, size=2, alpha=1)+
  scale_color_manual(values = friendly_pal("ito_seven")[c(3,6,7)], labels=c("Virginia", "Florida", "Africa"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(label=scientific)+
  labs(x="Chr. 3 Position", y=NULL, color="Pop.")+
  lims(y=c(-4,8))
  

  #facet_grid(pop~.)+
  #geom_vline(xintercept = 689841)
  #geom_hline(yintercept=5, linetype="dashed", color="grey50")
  


center<-plot_grid(fst.sc3, bp.sc3, ihs.sc3, rel_heights =c(0.3,0.3,0.4), ncol=1, labels=c("D", "E", "F"), align="v", axis="lr")

#lets also make plots of the region on scaffold_5 that has a peak for all 3 and is a CYP6a (insecticide resistance)

#ihs.VAall[CHR=="Scaffold_5"][order(IHS, decreasing=T)]
#peak is Scaffold_5  7875359

#make a plot centered on this

fst.sc5<-ggplot(z[chromosome=="Scaffold_5"& position<8600000&position>7500000])+geom_point(aes(x=position, y=fst.snp))+
  #scale_color_manual(values = friendly_pal("ito_seven")[3])+
  labs(x=NULL, y=NULL)+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  scale_y_continuous(limits=c(0,1), breaks=c(0, .25, .5, .75))


bp.sc5<-ggplot(bp.va.fl[CHR=="Scaffold_5"& pos<8600000&pos>7500000])+
  geom_point(aes(x=pos, y=M_XtX))+
  #scale_color_manual(values = friendly_pal("ito_seven")[3])+
  labs(x=NULL, y=NULL)+
  #theme(axis.text.x=element_blank())+
  guides(color="none")+
  theme(axis.text.x=element_blank())+
  lims(y=c(0,9))



ihs.sc5<-ggplot()+
  geom_point(data=ihs[scaffold=="Scaffold_5"&POSITION<8600000&POSITION>7500000&pop!="North America"], aes(x=POSITION, y=IHS, color=pop), alpha=0.5)+
  #geom_point(color="white", pch=21, size=2, alpha=1)+
  scale_color_manual(values = friendly_pal("ito_seven")[c(3,6,7)], labels=c("Virginia", "Florida", "Africa"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(label=scientific, breaks=c(7.5e6, 8e6,8.5e6))+
  labs(x="Chr. 5 Position", y=NULL, color="Pop.")+
  lims(y=c(-4,8))

right<-plot_grid(fst.sc5, bp.sc5, ihs.sc5, rel_heights =c(0.3,0.3,0.4), ncol=1, labels=c("G", "H", "I"), align="v", axis="lr")



jpeg("/scratch/perickso/private/ind_seq/Figures/Figure4_revision1.jpeg",  height=8, width=12, res=600, units="in")
plot_grid(left, center, right, nrow=1, align="v", axis="b")
dev.off()







