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
snptags.africa<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_Africa_snpstoplot.csv")


#FST
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsAfrica.Rdata")
fst.africa<-z



fst.plot.africa<-ggplot()+
  geom_rect(data=snptags.africa, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=fst.africa, aes(x=snp.id, y=fst.snp, color=chromosome))+
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

bp.plot.africa<-ggplot()+  
  geom_rect(data=snptags.africa, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=bp.va.af, aes(x=old.snp.id, y=M_XtX, color=as.factor(Scaffold)))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y="XtX")+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,15))

#IHS

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VAall<-as.data.table(wgscan.ihs$ihs)

setnames(info, "pos", "POSITION")
ihs.VAall<-merge(ihs.VAall, info, by=c("CHR", "POSITION"))
lab1<-ihs.VAall[CHR=="Scaffold_2"][order(IHS, decreasing=T)]$old.snp.id[1]
lab2<-ihs.VAall[CHR=="Scaffold_5"][order(IHS, decreasing=T)]$old.snp.id[1]

ihs.plot.africa<-ggplot()+
  geom_rect(data=snptags.africa, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=ihs.VAall,aes(x=old.snp.id, y=IHS, color=CHR))+
  scale_color_manual(values = friendly_pal("ito_seven"), labels=c("1", "2", "3", "4", "5"))+
  labs(x="SNP #", y="IHS", color="Chr.")+
  scale_x_continuous(label=scientific)+
  lims(y=c(-4,8))+
  guides(color="none")+
  annotate("text", x = lab1, y = 6, label = "*", color=friendly_pal("ito_seven")[2], size=12)+
  annotate("text", x = lab2, y = 6, label = "*", color=friendly_pal("ito_seven")[5], size=12)


  

left<-plot_grid(fst.plot.africa, bp.plot.africa, ihs.plot.africa, nrow=3, rel_heights=c(1, 1, 1.2), align="v", axis="lr", labels=c("a", "c", "e"))


###################
# now with FL data
###################

snptags.FL<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_snpstoplot_VA_ALL.csv")
load("/scratch/perickso/private/ind_seq/popgen/FST_VA_ALLvsFL.Rdata")
fst.FL<-z

fst.plot.FL<-ggplot()+
  geom_rect(data=snptags.FL, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=fst.FL, aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y=expression("Virginia-Florida F"[ST]))+
  theme(axis.text.x=element_blank())+  guides(color="none")+
  scale_y_continuous(limits=c(0,1), breaks=c(0, .25, .5, .75))

bp.va.fl<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/updated_data_for_manhattan.txt")
bp.va.fl[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]
setnames(info, "POSITION", "pos")
bp.va.fl<-merge(bp.va.fl, info, by=c("pos", "Scaffold"))

bp.plot.FL<-ggplot()+  
  geom_rect(data=snptags.FL, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=bp.va.fl, aes(x=old.snp.id, y=M_XtX, color=as.factor(Scaffold)))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y="XtX")+
  theme(axis.text.x=element_blank())+
  guides(color="none")+
  lims(y=c(0,15))



ihs.plot.FL<-ggplot()+
  geom_rect(data=snptags.FL, aes(xmin=first.snp, xmax=last.snp, ymin=-Inf, ymax=Inf), color="grey80")+
  geom_point(data=ihs.VAall,aes(x=old.snp.id, y=IHS, color=CHR))+
  scale_color_manual(values = friendly_pal("ito_seven"), labels=c("1", "2", "3", "4", "5"))+
  labs(x="SNP #", y="IHS", color="Chr.")+
  scale_x_continuous(label=scientific)+
  lims(y=c(-4,8))+
  guides(color="none")

right<-plot_grid(fst.plot.FL, bp.plot.FL, ihs.plot.FL, nrow=3, rel_heights=c(1, 1, 1.2), align="v", axis="lr", labels=c("b", "d", "f"))

#pdf("/scratch/perickso/private/ind_seq/Figures/Figure_5_Africa_FL_selection.pdf",  height=8, width=10)
#setEPS()
#postscript("/scratch/perickso/private/ind_seq/Figures/Figure_5_Africa_FL_selection.eps")

#jpeg("/scratch/perickso/private/ind_seq/Figures/Figure_5_Africa_FL_selection.jpeg",  height=8, width=10, res=600, units="in")
tiff("/scratch/perickso/private/ind_seq/Figures/Figure_5_Africa_FL_selection.tif",  height=8, width=10, res=600, units="in", compression="lzw")

plot_grid(left,right, nrow=1, align="h", axis="tb")
dev.off()

