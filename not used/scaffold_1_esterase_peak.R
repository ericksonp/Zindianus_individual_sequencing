



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




