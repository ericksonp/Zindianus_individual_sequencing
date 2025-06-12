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

bp.peaks<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_outlier_SNPs_q999.xtx")
setnames(bp.peaks, "SNP.info[outlier_SNPs_999$MRK, ]", "snp")
bp.peaks[,CHR:=paste0(tstrsplit(snp, split="_")[[1]],"_", tstrsplit(snp, split="_")[[2]]) ]
bp.peaks[,POSITION:=as.integer(tstrsplit(snp, split="_")[[3]])]
bp.peaks[order(M_XtX, decreasing=T)][CHR!="Scaffold_3"]
#Scaffold_1 20306622
#Scaffold_2 20748414
#Scaffold_4  6573623
#Scaffold_1 16187187
#Scaffold_4 16209741
#Scaffold_2_26766536
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

#what's within 10kb of each bayepass peak?

ann[chr=="Scaffold_1"&end>20306622-10000&start<20306622+10000&info=="gene"]
#adipocyte plasma protein
#unk (unkempt) -eye development

#this looks like it is close ot the scaffold 1 inversion breakpoint

ann[chr=="Scaffold_1"&end>16187187-10000&start<16187187+10000&info=="gene"]
#esterase-B1--insectiside resistance in Culex!!
# Scaffold_1	16177203 is a mis-sense variant in ANN01440 but is alanine to valine so not a big change
#maybe use vcfR to look for all NS variants in these genes

ann[chr=="Scaffold_2"&end>20748414-10000&start<20748414+10000&info=="gene"]

#amnionless - protein resporbtion and vitamin b12 in nephrocytes
#smoothened - hedgehog signaling

#this SNP is very close to the scaffold 2 inversion breakpoint

ann[chr=="Scaffold_2"&end>26766536-10000&start<26766536+10000&info=="gene"]
#no annotations

ann[chr=="Scaffold_4"&end>6573623-10000&start<6573623+10000&info=="gene"]
#rho guanine nucleotide exchange factor


ann[chr=="Scaffold_4"&end>16209741-10000&start<16209741+10000&info=="gene"]
#tbox transcription factor

#do the bayepass peaks come close to the inversion breakpoints?

inv<-fread("/scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed", header=F)


#what about the scaffold 5 location that is a peak for all 3 selection tests?
wins<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows.csv")
wins[CHR==5&all.3==T]
ann[chr=="Scaffold_5"&end>7844001-10000&start<7885001+10000&info=="gene"]

#this is a cluster of CYP450 genes

#read in snpeff data
eff<-fread("/scratch/perickso/private/ind_seq/SnpEff_annotations_data_table.csv")
eff[(gene=="ANN03537"|gene=="ANN03538"|gene=="ANN03539"|gene=="ANN03540")&type=="missense_variant"]
eff[(gene=="ANN03537")&type=="missense_variant"]

#IHS peak is scaffold 5 7875359
eff[POS==7875359] #synonymous variant but several missense variants nearby
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")
z[position>7844001&position<7886000&chromosome=="Scaffold_5"][order(fst.snp, decreasing=T)]
eff[POS==7844664]




## look at haplotypes

# hap_file="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.2.vcf" #change to ancestral after it's made
# hh <- data2haplohh(hap_file = hap_file,
#                    chr.name = "Scaffold_2",
#                    min_perc_geno.mrk = 90,
#                    polarize_vcf = TRUE,
#                    vcf_reader = "vcfR")
# samps<-data.table(id=hap.names(hh))
# samps[,haplotype.id:=1:(nrow(samps))]
# samps[,sample.id:=substr(id,1,nchar(id)-2 )]
# samps[,hapnum:=rep(c(1:2), times=.N/2)]
# metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
# 
# samps<-merge(metadata, samps, by="sample.id")
# haps.to.use<-samps[, haplotype.id] #can add other conditionals here for specific subsets
# 
# 
# hh_subset = subset(hh, select.hap = haps.to.use, min_perc_geno.mrk = 75, min_maf=0)
# 
# haplos<-as.data.table(t(hh_subset@haplo))
# haplos[,pos:=hh_subset@positions]
# 
# metadata[,index:=1:nrow(metadata)]
# haplos.melt<-melt(haplos, id.vars="pos",value.vars=names(haplos)[1:(ncol(haplos)-1)], variable.name="haplo.id", value.name="genotype" )
# haplos.melt[,sample.id:=gsub('.{2}$', '', haplo.id)]
# haplos.melt<-merge(haplos.melt, metadata[,.(sample.id, loc.spec)], by="sample.id")
# 
# save(haplos.melt, file="/scratch/perickso/private/ind_seq/popgen/scaffold_2haplotype_table.Rdat")
# 
# load("/scratch/perickso/private/ind_seq/popgen/scaffold_1haplotype_table.Rdat")
# 

##################################################
#Scaffold_5 plot
##################################################

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#load IHS data and choose SNPs
load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[CHR=="Scaffold_5"][order(IHS, decreasing=T)] # 7875359

ihs.high<-ihs.VA[CHR=="Scaffold_5"&POSITION>7840000&POSITION<7920000, POSITION]


scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
setnames(scaffolds, "V1" ,"chr")
scaffolds[,index:=1:nrow(scaffolds)]
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
# save(hh_subset, file="/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf5.Rdata")

load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf5.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf5.Rdata")

i=5
pos=7875359
mk<-scan.dt[POSITION==pos, markernum]
ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk)$ehh)
ehh.melt5_7<-melt(ehh, id.vars="POSITION")
save(ehh.melt5_7, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc5_7875359.Rdata")


ehhplot5<-ggplot(ehh.melt5_7)+
  geom_line(aes(x=POSITION/1000000, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 5 Position (Mb)", y="EHH", color=NULL)+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high)/1000000,max(ihs.high)/1000000) )+
  guides(color="none")


#make table that produces a relative position for each SNP
snp_spacing=round((max(ihs.high)-min(ihs.high))/length(ihs.high))
position.table<-data.table(pos=ihs.high,
                           count=c(1:(length(ihs.high))),
                           uniform_pos=seq(min(ihs.high), (max(ihs.high)-2*snp_spacing), by=snp_spacing),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
#can we add annotations to this?
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons<-ann[chr=="Scaffold_5"&info=="exon"&((start>min(ihs.high)&start<max(ihs.high))|(end>min(ihs.high)&end<max(ihs.high)))]
genes<-ann[chr=="Scaffold_5"&info=="gene"&((start>min(ihs.high)&start<max(ihs.high))|(end>min(ihs.high)&end<max(ihs.high)))]
#ANN03542 is not a cyp
#everything else is

top_scale<-ggplot()+
  geom_segment(data=genes, aes(x=start, xend=end, y=1.3, yend=1.3))+
  geom_segment(data=position.table, aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  #plot cyps in color
  geom_rect(data=exons[gene_name!="ANN03542"], aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5, fill=gene_name))+
  #plot other gene in gray
  geom_rect(data=exons[gene_name=="ANN03542"], aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5), fill="grey50")+
  scale_x_continuous(position="top", expand = c(0, 0), labels = label_scientific(base=10, digits=5), limits=c(min(ihs.high),max(ihs.high)) )+
  # labs(x="Chromosome 2 position")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  guides(fill="none")


 load("/scratch/perickso/private/ind_seq/popgen/scaffold_5haplotype_table.Rdat")

haplos.to.plot<-haplos.melt[pos%in%ihs.high]

focal.geno<-haplos.to.plot[pos==7875359] #ordering based on IHS peak
focal.geno[,fixed.geno:=genotype]
haplos.to.plot<-merge(haplos.to.plot, focal.geno[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot<-haplos.to.plot[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot[,haplo.index:=rleid(haplo.id)]
haplos.to.plot[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot<-unique(haplos.to.plot$pos)


haplo.plot<-ggplot(haplos.to.plot)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+

  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 



geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use<-info[chr=="Scaffold_5"&pos%in%snps.in.haplo.plot, snp.id]

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)

ggLD <- function(data){
  data <- tidyr::as_tibble(data)
  colnames(data) <- c(1:ncol(data))
  n <- length(data)
  
  # Tidy data, only taking unique pairs of data
  values <- data %>%
    dplyr::mutate(idx1 = c(1:nrow(data))) %>%
    tidyr::pivot_longer(!.data$idx1, names_to = "idx2", values_to = "LD") %>%
    dplyr::mutate(dplyr::across(idx2, as.double)) %>%
    dplyr::filter(!duplicated(paste(pmax(.data$idx1, .data$idx2), pmin(.data$idx1, .data$idx2), sep = "_"))) %>%
    tidyr::unite("id", .data$idx1:.data$idx2, remove = FALSE) %>%
    dplyr::mutate(diff = abs(idx2 - idx1))
  
  # Calculate coordinates for geom_polygon
  positions <- dplyr::bind_rows(values, values, values, values) %>%
    dplyr::group_by(diff, idx1) %>%
    dplyr::mutate(add_index1 = c(0, 1, 0, 1),
                  add_index2 = c(0, -1, 0, 1),
                  minus1_index = c(1, 1, 0, 1)) %>%
    dplyr::mutate(x = diff * 5 / n + 10 / n * (idx1 - minus1_index) + 5 / n * add_index1,
                  y = 5 - diff * 5 / n + 5 / n * add_index2) %>%
    dplyr::ungroup()
  
  # ggplot2
  positions %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_polygon(ggplot2::aes(fill = .data$LD, group = .data$id)) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_distiller(type = "seq", palette = 1, direction = 1)
}
ldmat<-ld$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt<-as.data.table(ldmat)
names(lddt)=as.character(ld$snp.id)
ld.plot<-ggLD(lddt)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 


jpeg("/scratch/perickso/private/ind_seq/Figures/scaffold_5_bpsignal_cpr.jpeg",  height=10, width=6, res=600, units="in")
plot_grid(ehhplot5, top_scale, haplo.plot, ld.plot, nrow=4, labels=c("A"," ", "B", "C" ), rel_heights=c(0.3,0.1, 0.4, 0.3), align="v",axis="lr")
dev.off()


#look for missense variants in highest IHS region

eff<-fread("/scratch/perickso/private/ind_seq/SnpEff_annotations_data_table.csv")
setnames(eff, c("CHROM", "POS"), c("CHR", "POSITION"))

eff<-merge(eff, ihs.VA, by=c("CHR", "POSITION"))
eff[CHR=="Scaffold_5"&IHS>4]
#IHS of 4.69 at 7883497, missesnse variant PHE->LEUC


##################################################
#Let's take a look at the IHS peak on scaffold 2
##################################################
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA[order(IHS, decreasing=T)][CHR=="Scaffold_2"]
#scaffold 2 26609601
# 
#  scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
#  setnames(scaffolds, "V1" ,"chr")
# scaffolds[,index:=1:nrow(scaffolds)]
# # 
# hap_file="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.2.vcf" #change to ancestral after it's made
# hh <- data2haplohh(hap_file = hap_file,
#                    chr.name = "Scaffold_2",
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
# save(hh_subset, file="/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf2.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf2.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf2.Rdata")

i=2
pos=26609601
mk<-scan.dt[POSITION==pos, markernum]
ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk)$ehh)
ehh.melt2_26<-melt(ehh, id.vars="POSITION")
save(ehh.melt2_26, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc2_26609601.Rdata")


load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")

ihs.high<-ihs.VA[CHR=="Scaffold_2"&POSITION>26590000&POSITION<26700000, POSITION]

ehhplot2<-ggplot(ehh.melt2_26)+
  geom_line(aes(x=POSITION/1000000, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 2 Position (Mb)", y="EHH", color=NULL)+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high)/1000000,max(ihs.high)/1000000) )+
  guides(color="none")

ann[chr=="Scaffold_2"&end>26609601-50000&start<26609601+100000&info=="gene"]
#this is cpr https://scijournals.onlinelibrary.wiley.com/doi/10.1002/ps.4852
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031037
#https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2583.2006.00647.x

load("/scratch/perickso/private/ind_seq/popgen/scaffold_2haplotype_table.Rdat")



#make table that produces a relative position for each SNP
snp_spacing=round((max(ihs.high)-min(ihs.high))/length(ihs.high))
position.table<-data.table(pos=ihs.high,
                           count=c(1:(length(ihs.high))),
                           uniform_pos=seq(min(ihs.high), (max(ihs.high)), by=snp_spacing),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
#can we add annotations to this?
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons<-ann[chr=="Scaffold_2"&info=="exon"&((start>min(ihs.high)&start<max(ihs.high))|(end>min(ihs.high)&end<max(ihs.high)))]
genes<-ann[chr=="Scaffold_2"&info=="gene"&((start>min(ihs.high)&start<max(ihs.high))|(end>min(ihs.high)&end<max(ihs.high)))]




top_scale<-ggplot()+
  geom_segment(data=genes, aes(x=start, xend=end, y=1.3, yend=1.3))+
  geom_segment(data=position.table, aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  geom_rect(data=exons, aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5, fill=ifelse(gene_name=="ANN06929", "black", "grey")))+
  scale_x_continuous(position="top", expand = c(0, 0), labels = label_scientific(base=10, digits=5), limits=c(min(ihs.high),max(ihs.high)) )+
 # labs(x="Chromosome 2 position")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c(friendly_pal("ito_seven")[2], "grey40"))+
  guides(fill="none")


#haplos.to.plot<-haplos.melt[pos%in%fst.high|pos%in%ihs.high]
haplos.to.plot<-haplos.melt[pos%in%ihs.high]

focal.geno<-haplos.to.plot[pos==26609601] #ordering based on IHS peak
focal.geno[,fixed.geno:=genotype]
haplos.to.plot<-merge(haplos.to.plot, focal.geno[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot<-haplos.to.plot[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot[,haplo.index:=rleid(haplo.id)]
haplos.to.plot[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot<-unique(haplos.to.plot$pos)


haplo.plot<-ggplot(haplos.to.plot)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+
  
  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 



geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use<-info[chr=="Scaffold_2"&pos%in%snps.in.haplo.plot, snp.id]

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)

ggLD <- function(data){
  data <- tidyr::as_tibble(data)
  colnames(data) <- c(1:ncol(data))
  n <- length(data)
  
  # Tidy data, only taking unique pairs of data
  values <- data %>%
    dplyr::mutate(idx1 = c(1:nrow(data))) %>%
    tidyr::pivot_longer(!.data$idx1, names_to = "idx2", values_to = "LD") %>%
    dplyr::mutate(dplyr::across(idx2, as.double)) %>%
    dplyr::filter(!duplicated(paste(pmax(.data$idx1, .data$idx2), pmin(.data$idx1, .data$idx2), sep = "_"))) %>%
    tidyr::unite("id", .data$idx1:.data$idx2, remove = FALSE) %>%
    dplyr::mutate(diff = abs(idx2 - idx1))
  
  # Calculate coordinates for geom_polygon
  positions <- dplyr::bind_rows(values, values, values, values) %>%
    dplyr::group_by(diff, idx1) %>%
    dplyr::mutate(add_index1 = c(0, 1, 0, 1),
                  add_index2 = c(0, -1, 0, 1),
                  minus1_index = c(1, 1, 0, 1)) %>%
    dplyr::mutate(x = diff * 5 / n + 10 / n * (idx1 - minus1_index) + 5 / n * add_index1,
                  y = 5 - diff * 5 / n + 5 / n * add_index2) %>%
    dplyr::ungroup()
  
  # ggplot2
  positions %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_polygon(ggplot2::aes(fill = .data$LD, group = .data$id)) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_distiller(type = "seq", palette = 1, direction = 1)
}
ldmat<-ld$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt<-as.data.table(ldmat)
names(lddt)=as.character(ld$snp.id)
ld.plot<-ggLD(lddt)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 


jpeg("/scratch/perickso/private/ind_seq/Figures/scaffold_2_bpsignal_cpr.jpeg",  height=10, width=6, res=600, units="in")
plot_grid(ehhplot2, top_scale, haplo.plot, ld.plot, nrow=4, labels=c("A"," ", "B", "C" ), rel_heights=c(0.3,0.1, 0.4, 0.3), align="v",axis="lr")
dev.off()









##################################################

#### now do this for scaffold 1 esterase peak####

##################################################

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

load("/scratch/perickso/private/ind_seq/popgen/scaffold_1haplotype_table.Rdat")
#Scaffold_1 16187187

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")

ihs.high<-ihs.VA[CHR=="Scaffold_1"&POSITION>16150000&POSITION<16250000, POSITION]

#make table that produces a relative position for each SNP
snp_spacing=round((max(ihs.high)-min(ihs.high))/length(ihs.high))
position.table<-data.table(pos=ihs.high,
                           count=c(1:(length(ihs.high))),
                           uniform_pos=seq(min(ihs.high), (max(ihs.high)+5*snp_spacing), by=snp_spacing),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
#can we add annotations to this?
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons<-ann[gene_name%in%c("ANN01440","ANN01441", "ANN01442", "ANN01443", "ANN01444") &info=="exon"]


top_scale<-ggplot(position.table)+
  geom_segment(aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  geom_rect(data=exons, aes(xmin=start,xmax=end, ymin=1.1, ymax=1.15, fill=gene_name))+
  scale_x_continuous(position="top", label=scientific, expand = c(0, 0))+
  labs(x="Chromosome 1 position")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())+
  scale_fill_manual(values = friendly_pal("bright_seven"))+
  guides(fill="none")


#haplos.to.plot<-haplos.melt[pos%in%fst.high|pos%in%ihs.high]
haplos.to.plot<-haplos.melt[pos%in%ihs.high]

focal.geno<-haplos.to.plot[pos==16187187] #ordering based on IHS peak
focal.geno[,fixed.geno:=genotype]
haplos.to.plot<-merge(haplos.to.plot, focal.geno[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot<-haplos.to.plot[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot[,haplo.index:=rleid(haplo.id)]
haplos.to.plot[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot<-unique(haplos.to.plot$pos)


haplo.plot<-ggplot(haplos.to.plot)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+
  
  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 



geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use<-info[chr=="Scaffold_1"&pos%in%snps.in.haplo.plot, snp.id]

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)

ggLD <- function(data){
  data <- tidyr::as_tibble(data)
  colnames(data) <- c(1:ncol(data))
  n <- length(data)
  
  # Tidy data, only taking unique pairs of data
  values <- data %>%
    dplyr::mutate(idx1 = c(1:nrow(data))) %>%
    tidyr::pivot_longer(!.data$idx1, names_to = "idx2", values_to = "LD") %>%
    dplyr::mutate(dplyr::across(idx2, as.double)) %>%
    dplyr::filter(!duplicated(paste(pmax(.data$idx1, .data$idx2), pmin(.data$idx1, .data$idx2), sep = "_"))) %>%
    tidyr::unite("id", .data$idx1:.data$idx2, remove = FALSE) %>%
    dplyr::mutate(diff = abs(idx2 - idx1))
  
  # Calculate coordinates for geom_polygon
  positions <- dplyr::bind_rows(values, values, values, values) %>%
    dplyr::group_by(diff, idx1) %>%
    dplyr::mutate(add_index1 = c(0, 1, 0, 1),
                  add_index2 = c(0, -1, 0, 1),
                  minus1_index = c(1, 1, 0, 1)) %>%
    dplyr::mutate(x = diff * 5 / n + 10 / n * (idx1 - minus1_index) + 5 / n * add_index1,
                  y = 5 - diff * 5 / n + 5 / n * add_index2) %>%
    dplyr::ungroup()
  
  # ggplot2
  positions %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_polygon(ggplot2::aes(fill = .data$LD, group = .data$id)) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_distiller(type = "seq", palette = 1, direction = 1)
}
ldmat<-ld$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt<-as.data.table(ldmat)
names(lddt)=as.character(ld$snp.id)
ld.plot<-ggLD(lddt)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) 


jpeg("/scratch/perickso/private/ind_seq/Figures/scaffold_1_bpsignal_esterase.jpeg",  height=10, width=6, res=600, units="in")
plot_grid(top_scale, haplo.plot, ld.plot, nrow=3, labels=c("", "A", "B" ), rel_heights=c(0.3, 0.4, 0.3), align="v",axis="lr")
dev.off()




