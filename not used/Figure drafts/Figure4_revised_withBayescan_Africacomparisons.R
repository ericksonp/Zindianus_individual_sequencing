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
snptags<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_Africa_snpstoplot.csv")


#FST
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsAfrica.Rdata")

fst.plot<-ggplot()+
  geom_rect(data=snptags, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=z, aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y=expression("Virginia-Africa F"[ST]))+
  theme(axis.text.x=element_blank())+  guides(color="none")+
  scale_y_continuous(limits=c(0,1), breaks=c(0, .25, .5, .75))



#Bayescan

bp.va.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/updated_data_for_manhattan.txt")
bp.va.af[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(old.snp.id=a$snp.id,
                 CHR=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,Scaffold:=as.integer(tstrsplit(CHR, split="_")[[2]])]
bp.va.af<-merge(bp.va.af, info, by=c("pos", "Scaffold"))

bp.plot<-ggplot()+  
  geom_rect(data=snptags, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=bp.va.af, aes(x=old.snp.id, y=M_XtX, color=as.factor(Scaffold)))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y="XtX")+
  theme(axis.text.x=element_blank())+
  guides(color="none")

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


jpeg("/scratch/perickso/private/ind_seq/Figures/Figure4_Africa_version.jpeg",  height=10, width=8, res=600, units="in")
left
dev.off()
