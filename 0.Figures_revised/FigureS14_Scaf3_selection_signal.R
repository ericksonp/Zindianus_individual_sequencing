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



##################################################
#Scaffold_3 plot
##################################################

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

#load IHS data and choose SNPs
load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[CHR=="Scaffold_3"][order(IHS, decreasing=T)] # 973443
min(ihs.VA[CHR=="Scaffold_3"&!is.na(IHS),POSITION])

ihs.high3<-ihs.VA[CHR=="Scaffold_3"&POSITION>400000&POSITION<1500000, POSITION]

bp.va.fl<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/updated_data_for_manhattan.txt")
bp.va.fl[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]
#get quantile from simulations
sims<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_sim_summary_pi_xtx.out")
q999<-quantile(sims$M_XtX,0.999)

bp.plot3<-ggplot()+geom_point(data=bp.va.fl[Scaffold==3&pos>400000&pos<1500000], aes(x=pos/1000000, y=M_XtX))+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high3)/1000000,max(ihs.high3)/1000000) )+
  theme(        axis.text.x=element_blank(),
                axis.title.x=element_blank())+
  labs(y="XtX")+
  geom_hline(aes(yintercept=q999), linetype="dashed", color="grey50")+
  lims(y=c(0,12.5))
  

load("/scratch/perickso/private/ind_seq/popgen/ehh_sc3_973443.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/CM_HPO_ehh_scan_scaf3.Rdata")
load("/scratch/perickso/private/ind_seq/popgen/hh_subset_scaf3.Rdata")


ehhplot3<-ggplot(ehh.melt1)+
  geom_line(aes(x=POSITION/1000000, y=value, color=variable))+
  scale_color_manual(values = friendly_pal("ito_seven")[c(4,6)],labels=c("Allele 1", "Allele 2") )+
  labs(x="Chr. 3 Position (Mb)", y="EHH", color=NULL)+
  scale_x_continuous(expand = c(0, 0),limits=c(min(ihs.high3)/1000000,max(ihs.high3)/1000000) )+
  theme(legend.position = c(0.65, 0.7))


#make table that produces a relative position for each SNP
snp_spacing3=round((max(ihs.high3)-min(ihs.high3))/length(ihs.high3))
position.table3<-data.table(pos=ihs.high3,
                           count=c(1:(length(ihs.high3))),
                           uniform_pos=seq(min(ihs.high3), (max(ihs.high3)-snp_spacing3), by=snp_spacing3),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
#can we add annotations to this?
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]
exons3<-ann[chr=="Scaffold_3"&info=="exon"&((start>min(ihs.high3)&start<max(ihs.high3))|(end>min(ihs.high3)&end<max(ihs.high3)))]
genes3<-ann[chr=="Scaffold_3"&info=="gene"&((start>min(ihs.high3)&start<max(ihs.high3))|(end>min(ihs.high3)&end<max(ihs.high3)))]

top_scale3<-ggplot()+
  geom_segment(data=genes3, aes(x=start, xend=end, y=1.3, yend=1.3))+
  geom_segment(data=position.table3, aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  #plot other gene in gray
  geom_rect(data=exons3, aes(xmin=start,xmax=end, ymin=1.1, ymax=1.5), fill="grey50")+
  scale_x_continuous(expand = c(0, 0), labels = label_scientific(base=10, digits=5), limits=c(min(ihs.high3),max(ihs.high3)) )+
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

load("/scratch/perickso/private/ind_seq/popgen/scaffold_3haplotype_table.Rdat")

haplos.to.plot3<-haplos.melt[pos%in%ihs.high3]

focal.geno3<-haplos.to.plot3[pos==973443] #ordering based on IHS peak
focal.geno3[,fixed.geno:=genotype]
haplos.to.plot3<-merge(haplos.to.plot3, focal.geno3[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot3<-haplos.to.plot3[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot3[,haplo.index:=rleid(haplo.id)]
haplos.to.plot3[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot3<-unique(haplos.to.plot3$pos)


haplo.plot3<-ggplot(haplos.to.plot3)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
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
samps.to.use3<-metadata[continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use3<-info[chr=="Scaffold_3"&pos%in%snps.in.haplo.plot3, snp.id]

ld3<-snpgdsLDMat(geno, sample.id=samps.to.use3, snp.id=snps.to.use3, method="composite", slide=0)

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
ldmat3<-ld3$LD^2 #need to square because "composite" returns correlation coefficient
#ldmat[upper.tri(ldmat)] <- NA
lddt3<-as.data.table(ldmat3)
names(lddt3)=as.character(ld3$snp.id)
ld.plot3<-ggLD(lddt3)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title.align=0.5) 




jpeg("/scratch/perickso/private/ind_seq/Figures/FigureS13_scaffold3_haplotypes.jpeg",  height=10, width=6, res=600, units="in")
plot_grid(bp.plot3, ehhplot3, top_scale3, haplo.plot3, ld.plot3, nrow=5, labels=c("a","b ", "", "c", "d" ), rel_heights=c(1.5,2,1,4,3), align="v",axis="lr")
dev.off()
