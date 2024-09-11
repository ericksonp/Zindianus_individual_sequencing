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



#FST
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")

fst.plot<-ggplot(z)+geom_point(aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x=NULL, y=expression("Florida-Virginia F"[ST]))+
  theme(axis.text.x=element_blank())+  guides(color="none")


load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VAall<-as.data.table(wgscan.ihs$ihs)

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(old.snp.id=a$snp.id,
                 CHR=a$chromosome,
                 POSITION=a$pos,
                 freq=a$afreq)
ihs.VAall<-merge(ihs.VAall, info, by=c("CHR", "POSITION"))


ihs.plot<-ggplot(ihs.VAall)+geom_point(aes(x=old.snp.id, y=IHS, color=CHR))+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x="\nGenome position", y="IHS")+
  theme(axis.text.x=element_blank())+
  guides(color="none")


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

fst.sc3<-ggplot(z[chromosome=="Scaffold_3"& position<2000000])+geom_point(aes(x=position, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("ito_seven")[3])+
  labs(x=NULL, y=expression("FL-VA F"[ST]))+
         theme(axis.text.x=element_blank())+
  guides(color="none")
  


ihs.sc3<-ggplot()+
  geom_point(data=ihs[scaffold=="Scaffold_3"&POSITION<2000000], aes(x=POSITION, y=IHS, color=pop), size=0.5)+
  scale_color_manual(values = friendly_pal("ito_seven")[c(7,6,5,4)])+
  theme(legend.position = "none")+
  labs(x="Chr. 3 position", y="IHS")+
  scale_x_continuous(label=scientific)+
  facet_grid(pop~.)+
  #geom_vline(xintercept = 689841)
  geom_hline(yintercept=5, linetype="dashed", color="grey50")

center<-plot_grid(fst.sc3, ihs.sc3, ncol=1, rel_heights=c(0.25, 0.75), labels=c("C", "D"), align="v", axis="lr")

#EHH

# scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
# setnames(scaffolds, "V1" ,"chr")
# scaffolds[,index:=1:nrow(scaffolds)]
# 
# hap_file="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.3.vcf" #change to ancestral after it's made
# hh <- data2haplohh(hap_file = hap_file,
#                    chr.name = "Scaffold_3",
#                    min_perc_geno.mrk = 90,
#                    polarize_vcf = TRUE,
#                    vcf_reader = "vcfR")
# samps<-data.table(id=hap.names(hh))
# samps[,haplotype.id:=1:(nrow(samps))]
# samps[,sample.id:=substr(id,1,nchar(id)-2 )]
# samps[,hapnum:=rep(c(1:2), times=.N/2)]
# 
# samps<-merge(metadata, samps, by="sample.id")
# haps.to.use<-samps[Location=="VA-CM"|Location=="VA-HPO", haplotype.id] #can add other conditionals here for specific subsets
# 
# hh_subset = subset(hh, select.hap = haps.to.use, min_perc_geno.mrk = 75, min_maf=0)
# 
# scan <- scan_hh(hh_subset)
# scan.dt<-as.data.table(scan)
# scan.dt[,markernum:=c(1:nrow(scan.dt))]


#IHS peak

i=3
pos=973443
# 
# mk<-scan.dt[POSITION==pos, markernum]
# ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
# ehh.melt1<-melt(ehh, id.vars="POSITION")
# save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")

load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")
ehhplot1<-ggplot(ehh.melt1)+
  geom_line(aes(x=POSITION, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 3 Position", y="EHH", color=NULL)+
  scale_x_continuous(label=scientific,  limits=c(400000,1250000), breaks=c(400000,800000, 1200000))+
  guides(color="none")

#FST peak

i=3
pos=689841
# load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")
# mk<-scan.dt[POSITION==pos, markernum]
# ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
# ehh.melt2<-melt(ehh, id.vars="POSITION")
# save(ehh.melt2, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc3_689841.Rdata")

load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_689841.Rdata")
ehhplot2<-ggplot(ehh.melt2)+
  geom_line(aes(x=POSITION, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x=NULL, y="EHH", color=NULL)+
  theme(axis.text.x=element_blank())+
  theme(legend.position=c(0.7, 0.85))+
  scale_x_continuous(limits=c(400000,1250000), breaks=c(400000,800000, 1200000))



left<-plot_grid(fst.plot, ihs.plot, nrow=2, align="v", axis="lr", labels=c("A", "B"))

right<-plot_grid(ehhplot2, ehhplot1, nrow=2, align="v", axis="lr", labels=c("E", "F"))


jpeg("/scratch/perickso/private/ind_seq/Figures/Figure4_revised.jpeg",  height=8, width=12, res=600, units="in")
plot_grid(left, center, right, nrow=1, align="v", axis="b")
dev.off()