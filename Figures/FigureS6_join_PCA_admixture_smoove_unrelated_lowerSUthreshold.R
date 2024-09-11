#at alan's suggestion, try recoloring PCA plots based on admixture plots--he noticed that a lot of the admixture individuals are split 50/50 between two groups, which could reflect inversions.


library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(lattice)
library(tidyr)
library(SeqArray)
library(stringr)
library(doMC)
registerDoMC(20)
library(lubridate)
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)

#all PCA
pca.by.chr<-fread("/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr.csv")

#North American PCA
pca.by.chr.na<-fread("/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_northamerica.csv")

#CM PCA
pca.by.chr.cm<-fread("/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_CMonly.csv")

#admixture data
y<-fread("/scratch/perickso/private/ind_seq/popgen/admixture_data_for_figure.csv")

# #assign k2 categories based on the proporition ancestry in group 2
# y[new.k.group==2&proportion<0.3, kprop:="<30%"]
# y[new.k.group==2&proportion>0.3&proportion<0.7, kprop:="30-70%"]
# y[new.k.group==2&proportion>0.7, kprop:=">70%"]
# 
# 
# y[new.k.group==3&proportion<0.3, kprop:="<30%"]
# y[new.k.group==3&proportion>0.3&proportion<0.7, kprop:="30-70%"]
# y[new.k.group==3&proportion>0.7, kprop:=">70%"]
# 
# y[new.k.group==4&proportion<0.3, kprop:="<30%"]
# y[new.k.group==4&proportion>0.3&proportion<0.7, kprop:="30-70%"]
# y[new.k.group==4&proportion>0.7, kprop:=">70%"]
# #combine PCA data and admixture data. new.k.group==2 is the data we are most interested in
# 
#  pca.by.chr.cm<-merge(pca.by.chr.cm, y[, .(new.k.group, sample.id, proportion, chr, kprop)], by=c("sample.id", "chr"), all=T)
#  pca.by.chr.na<-merge(pca.by.chr.na, y[, .(new.k.group,sample.id, proportion, chr, kprop)], by=c("sample.id", "chr"))
#  pca.by.chr<-merge(pca.by.chr, y[, .(new.k.group,sample.id, proportion, chr, kprop)], by=c("sample.id", "chr"))
# # 
# # 
#  
# # 
# # look at how the k group assignment affects pc1-pc2  
# ggplot(pca.by.chr.na[new.k.group==2])+geom_point(aes(x=pc1, y=pc2, color=kprop))+facet_wrap(~chr)
# #k2 groups on 1 and 2
# ggplot(pca.by.chr.na[new.k.group==3])+geom_point(aes(x=pc1, y=pc2, color=kprop))+facet_wrap(~chr)
# #groups on 1, 2, and 5
# ggplot(pca.by.chr.na[new.k.group==4])+geom_point(aes(x=pc1, y=pc2, color=kprop))+facet_wrap(~chr)
# #groups on 1


##Okay, so the admixture plots and PCs are totally correlated. Could these be the inversions previously found on scaffolds 1, 2, and 5? let's get those genotypes
library(SNPRelate)
library(vcfR)

vcf <- read.vcfR("/scratch/perickso/private/ind_seq/sv/zap_all_called_sv.smoove.square.vcf.gz")

#get data from INFO column with has info about the structural variantsinfo
sv.info<-as.data.table(INFO2df(vcf))

#extract fixed data into a data.table

v<-as.data.table(getFIX(vcf))

sv.data<-cbind(v, sv.info)
sv.data[,POS:=as.numeric(POS)]
sv.data[,SVLEN:=as.numeric(SVLEN)]

#assuming that any SVs large enough to influence PCs will be at least 500 kb
inv<-sv.data[as.numeric(SVLEN)>1000000&as.numeric(SU)>5, .(CHROM, POS)]

large<-sv.data[as.numeric(SVLEN)>100000&as.numeric(SU)>25, .(CHROM, POS)]
#9 on scaffold 1, 10 on scaffold 2, 1 on 3, 1 on 4, 7 on 5

setnames(inv, c("chr", "pos"))


#make data table of sv genotypes
geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/sv/zap_all_called_sv.smoove.square.vcf.gds" , allow.fork=T)

#make master table of infformation
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)

#set key for indexing
setkey(info, chr, pos)
#get snp.ids of the focal inversion
inv.ids<-unique(info[inv]$snp.id)

large.ids<-unique(info[large]$snp.id)

genos.large<-as.data.table(snpgdsGetGeno(geno, snp.id=large.ids))
genos.inv<-as.data.table(snpgdsGetGeno(geno, snp.id=inv.ids))


names(genos.inv)=paste0("inv", inv.ids)
names(genos.large)=paste0("sv", large.ids)

genos.inv[,sample.id:=snpgdsGetGeno(geno, snp.id=inv.ids, with.id=T)$sample.id]
genos.large[,sample.id:=snpgdsGetGeno(geno, snp.id=large.ids, with.id=T)$sample.id]


pca.by.chr.na.inv<-merge(pca.by.chr.na, genos.inv, by=c("sample.id"))
pca.by.chr.na.large<-merge(pca.by.chr.na, genos.large, by=c("sample.id"))

#make pdfs to loop through all svs and plot pcs

plot_data <- function(x){
  ggplot(pca.by.chr.na.inv, aes_string(x="pc1",y="pc2", color=x)) +
    geom_point() + facet_wrap(~chr)+
    ggtitle(paste(x))
}


pdf("/scratch/perickso/private/ind_seq/popgen/plots/test_pca_inv_newthreshold.pdf")
foreach(variant=names(genos.inv)[1:(ncol(genos.inv)-1)])%do%{
  plot_data(variant)
}
dev.off()


plot_data <- function(x){
  ggplot(pca.by.chr.na.large, aes_string(x="pc1",y="pc2", color=x)) +
    geom_point() + facet_wrap(~chr)+
    ggtitle(paste(x))
}


pdf("/scratch/perickso/private/ind_seq/popgen/plots/test_pca_sv.pdf")
foreach(variant=names(genos.large)[1:(ncol(genos.large)-1)])%do%{
  plot_data(variant)
}
dev.off()



ggplot(pca.by.chr.na.large[chr==1])+
  geom_point(aes(x=pc1, y=pc2, color=as.factor(interaction(sv2107, sv1523))))

#marker 2107 is Scaffold_1  9134187 #5163474 BP inversion
#marker 1523 is  Scaffold_1  7100638 # 9647970 BP "duplication" probably not right but something going on.

#setnames(pca.by.chr.na.large, c("sv2107", "sv1523"), c("inv_9.1_Mb", "dup_7.1_Mb"))
pca.by.chr.na.large[sv2107==0, inv_9.1_Mb:="1/1"]
pca.by.chr.na.large[sv2107==1, inv_9.1_Mb:="0/1"]
pca.by.chr.na.large[sv2107==2, inv_9.1_Mb:="0/0"]
pca.by.chr.na.large[sv1523==0, dup_7.1_Mb:="1/1"]
pca.by.chr.na.large[sv1523==1, dup_7.1_Mb:="0/1"]
pca.by.chr.na.large[sv1523==2, dup_7.1_Mb:="0/0"]



ggplot(pca.by.chr.na.large[chr==1])+
  geom_point(aes(x=pc1, y=pc2, color=as.factor(interaction(dup_7.1_Mb, inv_9.1_Mb, sep = " - "))))+
  labs(x="PC1", y="PC2", color="combined genotype\nduplication @ 7.1 Mb -\ninversion @ 9.1 Mb")+
  scale_color_manual(values = friendly_pal("ito_seven"))
  

jpeg('/scratch/perickso/private/ind_seq/Figures/scaffold_1_sv_pca.jpg', height=5, width=8, units="in", res=600)

ggplot(pca.by.chr.na.large[chr==1])+
  geom_point(aes(x=pc1, y=pc2, shape=inv_9.1_Mb, color=dup_7.1_Mb), size=2)+
  labs(x="PC1", y="PC2", shape="sv @ 9.1 Mb", color="sv @ 7.1 Mb")+
  scale_color_manual(values = friendly_pal("ito_seven")[3:5])

dev.off()

#do stats

#correlation between PCs and inversion at 9.1
summary(lm(pc1~sv2107, data=pca.by.chr.na.large))
summary(lm(pc2~sv2107, data=pca.by.chr.na.large))

#correlation between PCs and duplication at 7.1
summary(lm(pc1~sv1523, data=pca.by.chr.na.large))
summary(lm(pc2~sv1523, data=pca.by.chr.na.large))

#get sample names for individuals with each sv in order to make samplots

genos.large[sv1523==0|sv2107==1|sv1523==1, sample.id]


bam.files2<-paste0("/scratch/perickso/private/ind_seq/RGSM_final_bams/", genos.large[sv1523==0|sv2107==1|sv1523==1, sample.id], ".RG.bam")
fwrite(list(bam.files2), file="/scratch/perickso/private/ind_seq/sv/scaffold1_twoSVs.csv")
