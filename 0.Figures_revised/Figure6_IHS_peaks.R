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

ihs.high5<-ihs.VA[CHR=="Scaffold_5"&POSITION>7840000&POSITION<7920000, POSITION]

bp.va.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/updated_data_for_manhattan.txt")
bp.va.af[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]
#get quantile from simulations
sims<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/zap_ind_2023_sim_summary_pi_xtx.out")
q999<-quantile(sims$M_XtX,0.999)

bp.plot5<-ggplot()+geom_point(data=bp.va.af[Scaffold==5&pos>7840000&pos<7920000], aes(x=pos/1000000, y=M_XtX))+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high5)/1000000,max(ihs.high5)/1000000) )+
  theme(        axis.text.x=element_blank(),
                axis.title.x=element_blank())+
  labs(y="XtX")+
  geom_hline(aes(yintercept=q999), linetype="dashed", color="grey50")+
  lims(y=c(0,12.5))
  
#scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
#setnames(scaffolds, "V1" ,"chr")
#scaffolds[,index:=1:nrow(scaffolds)]
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

#i=5
#pos=7875359
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk)$ehh)
#ehh.melt5_7<-melt(ehh, id.vars="POSITION")
#save(ehh.melt5_7, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc5_7875359.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/ehh_sc5_7875359.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf5.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf5.Rdata")


ehhplot5<-ggplot(ehh.melt5_7)+
  geom_line(aes(x=POSITION/1000000, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 5 Position (Mb)", y="EHH", color=NULL)+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high5)/1000000,max(ihs.high5)/1000000) )+
  theme(legend.position = c(0.65, 0.7))


#make table that produces a relative position for each SNP
snp_spacing5=round((max(ihs.high5)-min(ihs.high5))/length(ihs.high5))
position.table5<-data.table(pos=ihs.high5,
                           count=c(1:(length(ihs.high5))),
                           uniform_pos=seq(min(ihs.high5), (max(ihs.high5)-2*snp_spacing5), by=snp_spacing5),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
#can we add annotations to this?
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons5<-ann[chr=="Scaffold_5"&info=="exon"&((start>min(ihs.high5)&start<max(ihs.high5))|(end>min(ihs.high5)&end<max(ihs.high5)))]
genes5<-ann[chr=="Scaffold_5"&info=="gene"&((start>min(ihs.high5)&start<max(ihs.high5))|(end>min(ihs.high5)&end<max(ihs.high5)))]
#ANN03542 is not a cyp
#everything else is

top_scale5<-ggplot()+
  geom_segment(data=genes5, aes(x=start, xend=end, y=1.3, yend=1.3))+
  geom_segment(data=position.table5, aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  #plot cyps in color
  geom_rect(data=exons5[gene_name!="ANN03542"], aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5, fill=gene_name))+
  #plot other gene in gray
  geom_rect(data=exons5[gene_name=="ANN03542"], aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5), fill="grey50")+
  scale_x_continuous(expand = c(0, 0), labels = label_scientific(base=10, digits=5), limits=c(min(ihs.high5),max(ihs.high5)) )+
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
  guides(fill="none")+
  lims(y=c(0,1.5))

load("/scratch/perickso/private/ind_seq/popgen/scaffold_5haplotype_table.Rdat")

haplos.to.plot5<-haplos.melt[pos%in%ihs.high5]

focal.geno5<-haplos.to.plot5[pos==7875359] #ordering based on IHS peak
focal.geno5[,fixed.geno:=genotype]
haplos.to.plot5<-merge(haplos.to.plot5, focal.geno5[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot5<-haplos.to.plot5[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]

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



geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use5<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use5<-info[chr=="Scaffold_5"&pos%in%snps.in.haplo.plot5, snp.id]

ld5<-snpgdsLDMat(geno, sample.id=samps.to.use5, snp.id=snps.to.use5, method="composite", slide=0)

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
ldmat5<-ld5$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt5<-as.data.table(ldmat5)
names(lddt5)=as.character(ld5$snp.id)
ld.plot5<-ggLD(lddt5)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title.align=0.5) 




#jpeg("/scratch/perickso/private/ind_seq/Figures/scaffold_5_bpsignal_cpr.jpeg",  height=10, width=6, res=600, units="in")
right<-plot_grid(bp.plot5, ehhplot5, top_scale5, haplo.plot5, ld.plot5, nrow=5, labels=c("b","d ", "", "f", "h" ), rel_heights=c(1.5,2,1,4,3), align="v",axis="lr")


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

#load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
#ihs.VA[order(IHS, decreasing=T)][CHR=="Scaffold_2"]
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


#i=2
#pos=26609601
#mk<-scan.dt[POSITION==pos, markernum]
#ehh <- as.data.table(calc_ehh(hh_subset, mrk=mk)$ehh)
#ehh.melt2_26<-melt(ehh, id.vars="POSITION")
#save(ehh.melt2_26, file="/scratch/perickso/private/ind_seq/popgen/ehh_sc2_26609601.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf2.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf2.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/ehh_sc2_26609601.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")

ihs.high2<-ihs.VA[CHR=="Scaffold_2"&POSITION>26590000&POSITION<26700000, POSITION]



bp.plot2<-ggplot()+geom_point(data=bp.va.af[Scaffold==2&pos>26590000&pos<26700000], aes(x=pos/1000000, y=M_XtX))+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high2)/1000000,max(ihs.high2)/1000000) )+
  theme(        axis.text.x=element_blank(),
                axis.title.x=element_blank())+
  labs(y="XtX")+
  geom_hline(aes(yintercept=q999), linetype="dashed", color="grey50")+
  lims(y=c(0,12.5))

ehhplot2<-ggplot(ehh.melt2_26)+
  geom_line(aes(x=POSITION/1000000, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 2 Position (Mb)", y="EHH", color=NULL)+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high2)/1000000,max(ihs.high2)/1000000) )+
  guides(color="none")

#ann[chr=="Scaffold_2"&end>26609601-50000&start<26609601+100000&info=="gene"]
#this is cpr https://scijournals.onlinelibrary.wiley.com/doi/10.1002/ps.4852
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031037
#https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2583.2006.00647.x




#make table that produces a relative position for each SNP
snp_spacing2=round((max(ihs.high2)-min(ihs.high2))/length(ihs.high2))
position.table2<-data.table(pos=ihs.high2,
                           count=c(1:(length(ihs.high2))),
                           uniform_pos=seq(min(ihs.high2), (max(ihs.high2)), by=snp_spacing2),
                           y1=0,
                           y2=1)



exons2<-ann[chr=="Scaffold_2"&info=="exon"&((start>min(ihs.high2)&start<max(ihs.high2))|(end>min(ihs.high2)&end<max(ihs.high2)))]
genes2<-ann[chr=="Scaffold_2"&info=="gene"&((start>min(ihs.high2)&start<max(ihs.high2))|(end>min(ihs.high2)&end<max(ihs.high2)))]


top_scale2<-ggplot()+
  geom_segment(data=genes2, aes(x=start, xend=end, y=1.3, yend=1.3))+
  geom_segment(data=position.table2, aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  geom_rect(data=exons2, aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5, fill=ifelse(gene_name=="ANN06929", "black", "grey")))+
  scale_x_continuous(expand = c(0, 0), labels = label_scientific(base=10, digits=5), limits=c(min(ihs.high2),max(ihs.high2)) )+
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
  guides(fill="none")+
  lims(y=c(0,1.5))

load("/scratch/perickso/private/ind_seq/popgen/scaffold_2haplotype_table.Rdat")


#haplos.to.plot<-haplos.melt[pos%in%fst.high|pos%in%ihs.high]
haplos.to.plot2<-haplos.melt[pos%in%ihs.high2]

focal.geno2<-haplos.to.plot2[pos==26609601] #ordering based on IHS peak
focal.geno2[,fixed.geno:=genotype]
haplos.to.plot2<-merge(haplos.to.plot2, focal.geno2[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot2<-haplos.to.plot2[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
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


samps.to.use2<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use2<-info[chr=="Scaffold_2"&pos%in%snps.in.haplo.plot2, snp.id]

ld2<-snpgdsLDMat(geno, sample.id=samps.to.use2, snp.id=snps.to.use2, method="composite", slide=0)

ldmat2<-ld2$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt2<-as.data.table(ldmat2)
names(lddt2)=as.character(ld2$snp.id)
ld.plot2<-ggLD(lddt2)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  guides(fill="none")


left<-plot_grid(bp.plot2, ehhplot2, top_scale2, haplo.plot2, ld.plot2, nrow=5, labels=c("a"," c","", "e", "g" ), rel_heights=c(1.5,2,1,4,3), align="v",axis="lr")


jpeg("/scratch/perickso/private/ind_seq/popgen/plots/selection_peaks_2and5.jpg", 
     height=10, 
     width=8,
     units="in", 
     res=300 )
plot_grid(left, right, nrow=1, align="h", axis="tb")
dev.off()




