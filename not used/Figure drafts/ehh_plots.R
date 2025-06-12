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

#decide on other SNPs to look at for EHH plots

wins<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows.csv")
wins[all.3==T]
wins[bp.fst==T]

#there were a few "all 3 windows on chromosome 5 @ ~7844001 - 7886000 bp
#get IHS data for this region

ihs.VAall[CHR=="Scaffold_5"&POSITION>7844001&POSITION<7886000][order(IHS)]
# 7875359 is IHS peak with 4.54 IHS

#use code below to get autosomal haplotype files for ehh

# scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
# setnames(scaffolds, "V1" ,"chr")
# scaffolds[,index:=1:nrow(scaffolds)]
# 
# hap_file="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.5.vcf" #change to ancestral after it's made
# hh <- data2haplohh(hap_file = hap_file,
#                    chr.name = "Scaffold_5",
#                    min_perc_geno.mrk = 90,
#                    polarize_vcf = TRUE,
#                    vcf_reader = "vcfR")
# samps<-data.table(id=hap.names(hh))
# samps[,haplotype.id:=1:(nrow(samps))]
# samps[,sample.id:=substr(id,1,nchar(id)-2 )]
# samps[,hapnum:=rep(c(1:2), times=.N/2)]
# metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
# 
# samps<-merge(metadata, samps, by="sample.id")
# haps.to.use<-samps[Location=="VA-CM"|Location=="VA-HPO", haplotype.id] #can add other conditionals here for specific subsets
# 
# hh_subset = subset(hh, select.hap = haps.to.use, min_perc_geno.mrk = 75, min_maf=0)
# 
# scan <- scan_hh(hh_subset)
# scan.dt<-as.data.table(scan)
# scan.dt[,markernum:=c(1:nrow(scan.dt))]
# save(scan.dt, file="/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf5.Rdata")

#this code is to get X chromosome haplotype files (adjusting males) for ehh
scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
setnames(scaffolds, "V1" ,"chr")
scaffolds[,index:=1:nrow(scaffolds)]

hap_file="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.3.vcf" #change to ancestral after it's made
hh <- data2haplohh(hap_file = hap_file,
                   chr.name = "Scaffold_3",
                   min_perc_geno.mrk = 90,
                   polarize_vcf = TRUE,
                   vcf_reader = "vcfR")
samps<-data.table(id=hap.names(hh))
samps[,haplotype.id:=1:(nrow(samps))]
samps[,sample.id:=substr(id,1,nchar(id)-2 )]
samps[,hapnum:=rep(c(1:2), times=.N/2)]
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)

samps<-merge(metadata, samps, by="sample.id")
#next line only grabs one haplotype for males
haps.to.use<-samps[(assigned_sex=="F"|(assigned_sex=="M"&hapnum==1))&(Location=="VA-CM"|Location=="VA-HPO"), haplotype.id] #can add other conditionals here for specific subsets 

hh_subset = subset(hh, select.hap = haps.to.use, min_perc_geno.mrk = 75, min_maf=0)

scan <- scan_hh(hh_subset)
scan.dt<-as.data.table(scan)
scan.dt[,markernum:=c(1:nrow(scan.dt))]
save(scan.dt, file="/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf3.Rdata")

#IHS peak scaf 5

#load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf5.Rdata")
#i=5
#pos=7875359
# 
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
#ehh.melt1<-melt(ehh, id.vars="POSITION")
#save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc5_7875359.Rdata")

#what about Bayepass peaks on autosomes that pass simulation cutoff?
bp.peaks<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_outlier_SNPs_q999.xtx")
setnames(bp.peaks, "SNP.info[outlier_SNPs_999$MRK, ]", "snp")
bp.peaks[,CHR:=paste0(tstrsplit(snp, split="_")[[1]],"_", tstrsplit(snp, split="_")[[2]]) ]
bp.peaks[,POSITION:=as.integer(tstrsplit(snp, split="_")[[3]])]

bp.peaks<-merge(ihs.VAall, bp.peaks, by=c("CHR", "POSITION"))
bp.peaks[order(M_XtX)]
#list of SNPs where we have IHS data and they are an outlier in XtX analysis
#Scaffold_1 16187187 --> this is the esterase peak
#Scaffold_2 26802966 
#Scaffold_4 16627234
#Scaffold_2 20593467 

#load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf2.Rdata")

#i=2
#pos=26802966
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
#ehh.melt1<-melt(ehh, id.vars="POSITION")
#save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc2_26802966.Rdata")


#i=2
#pos=20593467
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
#ehh.melt1<-melt(ehh, id.vars="POSITION")
#save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc2_20593467.Rdata")


#load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf1.Rdata")

#i=1
#pos=16187187
# 
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
#ehh.melt1<-melt(ehh, id.vars="POSITION")
#save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc1_16187187.Rdata")

#IHS peak scaf 3
load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf3.Rdata")

i=3
pos=973443

 mk<-scan.dt[POSITION==pos, markernum]
 ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
 ehh.melt1<-melt(ehh, id.vars="POSITION")
 save(ehh.melt1, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")
 
 #FST peak
 i=3
 pos=689841
  mk<-scan.dt[POSITION==pos, markernum]
  ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk )$ehh)
  ehh.melt2<-melt(ehh, id.vars="POSITION")
  save(ehh.melt2, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc3_689841.Rdata")
 
  
  
  ### Make plots from all of these data
 
load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")
ehhplot1<-ggplot(ehh.melt1)+
  geom_line(aes(x=POSITION, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 3 Position", y="EHH", color=NULL)+
  scale_x_continuous(label=scientific,  limits=c(400000,1250000), breaks=c(400000,800000, 1200000))+
  guides(color="none")


load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_689841.Rdata")
ehhplot2<-ggplot(ehh.melt2)+
  geom_line(aes(x=POSITION, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x=NULL, y="EHH", color=NULL)+
  theme(axis.text.x=element_blank())+
  theme(legend.position=c(0.7, 0.85))+
  scale_x_continuous(limits=c(400000,1250000), breaks=c(400000,800000, 1200000))