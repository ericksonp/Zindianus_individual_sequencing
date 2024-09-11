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
library(ggpubfigs)
library(dplyr)


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

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=tstrsplit(a$chromosome, split="_")[[2]],
                 pos=a$pos,
                 freq=a$afreq)
setkey(info, chr, pos)


metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[continent=="NorthAmerica", sample.id]

ff <- snpgdsSNPRateFreq(geno, sample.id=samps.to.use)
info[,missing:=ff$MissingRate]
info[,localfreq:=ff$AlleleFreq]

#set seed before random sampling
set.seed(43920238)

###############################
# Scaffold 1 #################
###############################
snps.to.use<-info[chr==1&missing==0&localfreq>0&localfreq<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2
lddt1<-as.data.table(ldmat)
names(lddt1)=as.character(ld$snp.id)


a<-ggLD(lddt1)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 1        North America")



###############################
# Scaffold 2 #################
###############################
snps.to.use<-info[chr==2&missing==0&localfreq>0&localfreq<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2
lddt2<-as.data.table(ldmat)
names(lddt2)=as.character(ld$snp.id)


b<-ggLD(lddt2)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma" , direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 2")+
  guides(fill="none")

###############################
# Scaffold 3 North America ####
###############################

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use.3<-metadata[continent=="NorthAmerica"&assigned_sex=="F", sample.id]

ff3 <- snpgdsSNPRateFreq(geno, sample.id=samps.to.use.3)
info[,missing.3:=ff3$MissingRate]
info[,localfreq.3:=ff3$AlleleFreq]



snps.to.use<-info[chr==3&missing.3==0&localfreq.3>0&localfreq.3<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use.3, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt3<-as.data.table(ldmat)
names(lddt3)=as.character(ld$snp.id)


c<-ggLD(lddt3)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 3")+
  guides(fill="none")



###############################
# Scaffold 4  #################
###############################
snps.to.use<-info[chr==4&missing==0&(localfreq>0&localfreq<1), snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt4<-as.data.table(ldmat)
names(lddt4)=as.character(ld$snp.id)


d<-ggLD(lddt4)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 4")+
  guides(fill="none")


###############################
# Scaffold 5 #################
###############################
snps.to.use<-info[chr==5&missing==0&(localfreq>0&localfreq<1), snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt5<-as.data.table(ldmat)
names(lddt5)=as.character(ld$snp.id)


e<-ggLD(lddt5)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma" , direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 5")+
  guides(fill="none")

########################
#Now repeat for Africa #
########################

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use<-metadata[continent=="Africa", sample.id]

ff <- snpgdsSNPRateFreq(geno, sample.id=samps.to.use)
info[,missing:=ff$MissingRate]
info[,localfreq:=ff$AlleleFreq]

#set seed before random sampling
set.seed(43920238)

###############################
# Scaffold 1 #################
###############################
snps.to.use<-info[chr==1&missing==0&localfreq>0&localfreq<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2
lddt6<-as.data.table(ldmat)
names(lddt6)=as.character(ld$snp.id)


aa<-ggLD(lddt6)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 1        Africa  ")



###############################
# Scaffold 2 #################
###############################
snps.to.use<-info[chr==2&missing==0&localfreq>0&localfreq<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2
lddt7<-as.data.table(ldmat)
names(lddt7)=as.character(ld$snp.id)


ba<-ggLD(lddt7)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 2")+
  guides(fill="none")

###############################
# Scaffold 3 North America ####
###############################

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
samps.to.use.3<-metadata[continent=="Africa"&assigned_sex=="F", sample.id]

ff3 <- snpgdsSNPRateFreq(geno, sample.id=samps.to.use.3)
info[,missing.3:=ff3$MissingRate]
info[,localfreq.3:=ff3$AlleleFreq]



snps.to.use<-info[chr==3&missing.3==0&localfreq.3>0&localfreq.3<1, snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use.3, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt8<-as.data.table(ldmat)
names(lddt8)=as.character(ld$snp.id)


ca<-ggLD(lddt8)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 3")+
  guides(fill="none")



###############################
# Scaffold 4  #################
###############################
snps.to.use<-info[chr==4&missing==0&(localfreq>0&localfreq<1), snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt9<-as.data.table(ldmat)
names(lddt9)=as.character(ld$snp.id)


da<-ggLD(lddt9)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 4")+
  guides(fill="none")


###############################
# Scaffold 5 #################
###############################
snps.to.use<-info[chr==5&missing==0&(localfreq>0&localfreq<1), snp.id]
snps.to.use<-sample(snps.to.use, 4000)
snps.to.use<-sort(snps.to.use)

ld<-snpgdsLDMat(geno, sample.id=samps.to.use, snp.id=snps.to.use, method="composite", slide=-1)

ldmat<-ld$LD^2

lddt10<-as.data.table(ldmat)
names(lddt10)=as.character(ld$snp.id)


ea<-ggLD(lddt10)+labs(fill="LD")+
  #theme(plot.margin = unit'(c(0,0.8,0,0.5), "cm"))+
  scale_fill_viridis(option="magma", direction = -1)+ 
  theme(legend.position = c(0.9, 0.4))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0))+
  labs(title="       Chr. 5")+
  guides(fill="none")



jpeg("/scratch/perickso/private/ind_seq/Figures/all_scaffold_inversion_NA_Africa.jpg", height=12, width=6, units="in", res=300 )

plot_grid(a,aa,b,ba,c,ca,d,da,e,ea, nrow=5,labels=c("A","A'", "B", "B'", "C", "C'", "D", "D'", "E", "E'"))

dev.off()
