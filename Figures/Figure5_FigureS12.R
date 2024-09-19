library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(viridis)
library(rehh)
library(ggsci)
library(ggLD)
library(ggpubfigs)
library(dplyr)
library(scales)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)

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

scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
setnames(scaffolds, "V1" ,"chr")
scaffolds[,index:=1:nrow(scaffolds)]


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
# haps.to.use<-samps[(assigned_sex=="F"|(assigned_sex=="M"&hapnum==1)), haplotype.id] #can add other conditionals here for specific subsets 
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
# save(haplos.melt, file="/scratch/perickso/private/ind_seq/popgen/scaffold_3haplotype_table.Rdat")
load("/scratch/perickso/private/ind_seq/popgen/scaffold_3haplotype_table.Rdat")




#plot for paper, using FST and IHS to choose which snps to plot

#load IHS data and choose SNPs
load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[order(IHS)] #973443

ihs.high<-ihs.VA[IHS>5&CHR=="Scaffold_3"&POSITION<1500000, POSITION]

#make table that produces a relative position for each SNP
snp_spacing=round((max(ihs.high)-min(ihs.high))/400)
position.table<-data.table(pos=ihs.high,
                           count=c(1:length(ihs.high)),
                           uniform_pos=seq(min(ihs.high), max(ihs.high), by=snp_spacing),
                           y1=0,
                           y2=1)


#make plot that shows spacing of SNPs
top_scale<-ggplot(position.table)+
  geom_segment(aes(x=uniform_pos, xend=pos, y=y1, yend=y2), linewidth=0.1)+
  scale_x_continuous(position="top", label=scientific, expand = c(0, 0))+
  labs(x="Chromosome 3 position")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())

#load FST

load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")


#load allele frequency differences

#haplos.to.plot<-haplos.melt[pos%in%fst.high|pos%in%ihs.high]
haplos.to.plot<-haplos.melt[pos%in%ihs.high]

focal.geno<-haplos.to.plot[pos==973443] #ordering based on IHS peak
focal.geno[,fixed.geno:=genotype]
haplos.to.plot<-merge(haplos.to.plot, focal.geno[,.(haplo.id, fixed.geno)], by="haplo.id")
haplos.to.plot<-haplos.to.plot[loc.spec%in%c("Africa", "VA-CM", "VA-HPO", "FL")][order(loc.spec, -fixed.geno)]
haplos.to.plot[,haplo.index:=rleid(haplo.id)]
haplos.to.plot[order(pos),pos.id:=rleid(pos)]

snps.in.haplo.plot<-unique(haplos.to.plot$pos)

#break in snpps from 632k to 840k which is position 11 in list, need to illustrate that in plot
  
haplo.plot<-ggplot(haplos.to.plot)+geom_tile(aes(x=pos.id, y=haplo.index, fill=as.factor(genotype)))+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_grid(loc.spec~.,scales="free_y", space="free_y", switch="y")+guides(fill="none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = friendly_pal("bright_seven")[c(5,2)])+
  #geom_vline(xintercept=64)+
  #geom_vline(xintercept=111)+
  theme(axis.line = element_blank())+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  geom_vline(xintercept=111)
  
  



#make LD heatmap for these snps:
geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[assigned_sex=="F"&continent=="NorthAmerica"&loc.spec!="FL", sample.id]


snps.to.use<-info[chr=="Scaffold_3"&pos%in%snps.in.haplo.plot, snp.id]

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)

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

#plot FST and IHS above

ihs.to.plot<-ihs.VA[POSITION%in%snps.in.haplo.plot&CHR=="Scaffold_3"]
ihs.to.plot[,index:=c(1:nrow(ihs.to.plot))]

ihs.plot<-ggplot()+geom_point(data=ihs.to.plot, 
                   aes(x=index, y=IHS,group=CHR))+
  labs(y="IHS", x=NULL)+
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())+
  lims(y=c(0,8))+
  scale_x_continuous(limits=c(0,400), expand = c(0, 0))+
  geom_vline(xintercept=111, linetype="dashed", color="gray50")


fst.to.plot<-z[position%in%snps.in.haplo.plot&chromosome=="Scaffold_3"]
fst.to.plot[,index:=c(1:nrow(fst.to.plot))]
fst.plot<-ggplot()+geom_point(data=fst.to.plot, 
                   aes(x=index, y=fst.snp, group=chromosome))+
  labs(y="FST", x=NULL)+
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())+
  geom_vline(xintercept=111, linetype="dashed", color="gray50")+
  scale_x_continuous(limits=c(0,400), expand = c(0, 0))


jpeg("/scratch/perickso/private/ind_seq/Figures/haplotype_plusLD_revised.jpeg", height=10, width=6,units="in", res=300 )

plot_grid(fst.plot, ihs.plot, haplo.plot, ld.plot, nrow=4, labels=c("A", "B", "C", "D"), rel_heights=c(0.1, 0.1, 0.5, 0.3), align="v",axis="lr")
dev.off()

jpeg("/scratch/perickso/private/ind_seq/Figures/haplotype_plusLD_revised_withscale.jpeg", height=10, width=6,units="in", res=300 )

plot_grid(top_scale, fst.plot, ihs.plot, haplo.plot, ld.plot, nrow=5, labels=c("", "A", "B", "C", "D"), rel_heights=c(0.1, 0.1, 0.1, 0.4, 0.3), align="v",axis="lr")
dev.off()


#plot LD in Africa and in FL--> supplemental figure

#africa
samps.to.use<-metadata[assigned_sex=="F"&continent=="Africa", sample.id]
snps.to.use<-info[chr=="Scaffold_3"&pos%in%snps.in.haplo.plot, snp.id]
ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)
ldmat<-ld$LD^2
#ldmat[upper.tri(ldmat)] <- NA
lddt<-as.data.table(ldmat)
names(lddt)=as.character(ld$snp.id)
af.ld.plot<-ggLD(lddt)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))


samps.to.use<-metadata[assigned_sex=="F"&loc.spec=="FL", sample.id]
snps.to.use<-info[chr=="Scaffold_3"&pos%in%snps.in.haplo.plot, snp.id]
ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=0)
ldmat<-ld$LD^2
#ldmat[upper.tri(ldmat)] <- NA
lddt<-as.data.table(ldmat)
names(lddt)=as.character(ld$snp.id)
fl.ld.plot<-ggLD(lddt)+labs(fill="LD")+
  #theme(plot.margin = unit(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))

jpeg("/scratch/perickso/private/ind_seq/popgen/plots/africa_florida_LD.jpg", height=5, width=4.5,units="in", res=300 )

plot_grid(af.ld.plot, fl.ld.plot, nrow=2, labels=c("A", "B"))
dev.off()



