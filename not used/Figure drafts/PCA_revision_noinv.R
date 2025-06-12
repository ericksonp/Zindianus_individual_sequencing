

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
library(ggsci)
library(ggpubfigs)

#make GDS
# snpgdsVCF2GDS(vcf.fn = "/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz",  
#             out.fn="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz.gds", 
#            method="biallelic.only")


#load genofile
genofile <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz.gds" , allow.fork=T)
samps <- read.gdsn(index.gdsn(genofile, "sample.id")) 

#get metadata
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
metadata[continent=="NorthAmerica", continent:="N. America"]
metadata[continent=="SouthAmerica", continent:="S. America"]
#exclude Kenya 2018 outlier which is way off on its own
samps<-samps[samps!="SRR14982088"]
#get snp info
a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]


#read in unrelated individuals
#unrelated<-fread("/scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id")
#names(unrelated)="id"

#get informative SNPs excluding X chromosome and filtering at MAF <0.05
snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          sample.id=samps,
                          snp.id=info[chr!="Scaffold_3", snp.id],
                          autosome.only=FALSE)


#turn SNPsets into vectors
snpset.id <-unlist(snpset)


 #calculate PCs

     pca.temp <- snpgdsPCA(genofile 
                           ,snp.id=snpset.id, 
                           autosome.only=FALSE, 
                           num.thread=20, 
                           sample.id=samps)

     #snpgdsPCASNPLoading(pca.temp, genofile)

     pca.dt <- data.table(pc1=pca.temp$eigenvect[,1],
                          pc2=pca.temp$eigenvect[,2],
                          pc3=pca.temp$eigenvect[,3],
                          pc4=pca.temp$eigenvect[,4],
                          sample.id=pca.temp$sample)
    


 pca.dt=merge(pca.dt, metadata, by="sample.id")
fwrite(pca.dt, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_mac3_noinv.csv")
 #plot overall PCA
 pca.percent <- data.table(pc1=pca.temp$eigenval[1]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc2=pca.temp$eigenval[2]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc3=pca.temp$eigenval[3]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc4=pca.temp$eigenval[4]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc5=pca.temp$eigenval[5]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc6=pca.temp$eigenval[6]/sum(pca.temp$eigenval, na.rm=T)*100)
                   
         

#calculate stats for North America-Africa PC1 divergence
 t.test(pca.dt[continent%in%c("N. America", "Africa")]$pc1~pca.dt[continent%in%c("N. America", "Africa")]$continent)
                

################################
#PCA only North American samples
###############################

#north american samples, autosomes only
#recalcuate ld pruning
snpset.na <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. America", sample.id]],
                             snp.id=info[chr!="Scaffold_3",snp.id],
                             maf=(3/(2*length(samps[samps%in%metadata[continent=="N. America", sample.id]]))))


snpset.id.na <-unlist(snpset.na)

pca.temp2 <- snpgdsPCA(genofile ,
                       snp.id=snpset.id.na, 
                       autosome.only=FALSE, 
                       sample.id=samps[samps%in%metadata[continent=="N. America", sample.id]],
                       num.thread=10)

pca.dt2 <- data.table(pc1=pca.temp2$eigenvect[,1],
                      pc2=pca.temp2$eigenvect[,2],
                      pc3=pca.temp2$eigenvect[,3],
                      pc4=pca.temp2$eigenvect[,4],
                      sample.id=pca.temp2$sample)
pca.percent2 <- data.table(pc1=pca.temp2$eigenval[1]/sum(pca.temp2$eigenval, na.rm=T)*100,
                           pc2=pca.temp2$eigenval[2]/sum(pca.temp2$eigenval, na.rm=T)*100,
                           pc3=pca.temp2$eigenval[3]/sum(pca.temp2$eigenval, na.rm=T)*100,
                           pc4=pca.temp2$eigenval[4]/sum(pca.temp2$eigenval, na.rm=T)*100)


pca.dt2=merge(pca.dt2, metadata, by="sample.id")
pca.dt2[,loc.spec:=factor(loc.spec, levels=c("Northeast", "VA-CM", "VA-HPO", "TN", "NC", "FL"))]
pca.dt2<-pca.dt2[order(loc.spec)]

fwrite(pca.dt2, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_northamerica_mac3_noinv.csv")



summary(aov(pca.dt2$pc1~pca.dt2$loc.spec))
summary(aov(pca.dt2$pc2~pca.dt2$loc.spec))


########################################
#Carter Mountain by year
########################################


snpset.cm <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=3/(2*length(samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]])))


snpset.id.cm <-unlist(snpset.cm)

pca.temp3 <- snpgdsPCA(genofile ,
                       snp.id=snpset.id.cm, 
                       autosome.only=FALSE, 
                       sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]],
                       num.thread=10)

pca.dt3 <- data.table(pc1=pca.temp3$eigenvect[,1],
                      pc2=pca.temp3$eigenvect[,2],
                      pc3=pca.temp3$eigenvect[,3],
                      pc4=pca.temp3$eigenvect[,4],
                      sample.id=pca.temp3$sample)

pca.percent3 <- data.table(pc1=pca.temp3$eigenval[1]/sum(pca.temp3$eigenval, na.rm=T)*100,
                           pc2=pca.temp3$eigenval[2]/sum(pca.temp3$eigenval, na.rm=T)*100,
                           pc3=pca.temp3$eigenval[3]/sum(pca.temp3$eigenval, na.rm=T)*100,
                           pc3=pca.temp3$eigenval[4]/sum(pca.temp3$eigenval, na.rm=T)*100)

metadata[Season=="mid", Season:="early"]
pca.dt3=merge(pca.dt3, metadata, by="sample.id")


fwrite(pca.dt3, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_CMonly_mac3_noinv.csv")



#stats
summary(aov(pca.dt3$pc1~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc2~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc3~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc4~as.factor(pca.dt3$Year)))


pca.dt[,continent:=factor(continent, levels=c("Africa", "N. America", "S. America", "Hawaii"))]   
pca.dt<-pca.dt[order(continent)]


pa<-ggplot(pca.dt, aes(x=pc1, y=pc2, fill=continent))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC1,", round(pca.percent[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent[1,2], 2), "%", sep=" "),
       fill="Continent")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  #theme(legend.position = "bottom")+
  guides(fill="none")

pap<-ggplot(pca.dt, aes(x=pc3, y=pc4, fill=continent))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC3,", round(pca.percent[1,3], 2), "%", sep=" "),
       y=paste("PC4,", round(pca.percent[1,4], 2), "%", sep=" "),
       fill="Continent")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2))

pb<-ggplot(pca.dt2, aes(x=pc1, y=pc2, fill=loc.spec))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC1,", round(pca.percent2[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent2[1,2], 2), "%", sep=" "),
       fill="Location")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  #theme(legend.position = "bottom")+
  guides(fill="none")


pbp<-ggplot(pca.dt2, aes(x=pc3, y=pc4, fill=loc.spec))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC3,", round(pca.percent2[1,3], 2), "%", sep=" "),
       y=paste("PC4,", round(pca.percent2[1,4], 2), "%", sep=" "),
       fill="Location")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  theme(legend.position = "bottom")

pc<-ggplot(pca.dt3, aes(x=pc1, y=pc2, fill=as.factor(Year)))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC1,", round(pca.percent3[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent3[1,2], 2), "%", sep=" "),
       fill="Year")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  #theme(legend.position = "bottom")+
  guides(fill="none")

pcp<-ggplot(pca.dt3, aes(x=pc3, y=pc4, fill=as.factor(Year)))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x=paste("PC3,", round(pca.percent3[1,3], 2), "%", sep=" "),
       y=paste("PC4,", round(pca.percent3[1,4], 2), "%", sep=" "),
       fill="Year")+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2))

#pdf("/scratch/perickso/private/ind_seq/Figures/Figure_1_PCA_mac3_noinv.pdf", height=4, width=12)
#plot_grid(pa, pb, pc, nrow=1, labels=c("A", "B", "C"))
#dev.off()

pdf("/scratch/perickso/private/ind_seq/Figures/Figure_1_PCA_mac3_noinv_pc34.pdf", height=8, width=12)
plot_grid(pa, pb, pc, pap,pbp,pcp, nrow=2, labels=c("A", "B", "C", "A'", "B'", "C'"), align="vh", axis="tblr")
dev.off()



##########################
### FOR SUPPLEMENT: repeat all analyses by chromosome
#do PCA by individual chromosomes; only use females for scaffold 3
#########################
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
library(ggsci)
library(ggpubfigs)




#load genofile
genofile <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz.gds" , allow.fork=T)
samps <- read.gdsn(index.gdsn(genofile, "sample.id")) 

#get metadata
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
metadata[continent=="NorthAmerica", continent:="N. America"]
metadata[continent=="SouthAmerica", continent:="S. America"]
#exclude Kenya 2018 outlier
samps<-samps[samps!="SRR14982088"]
#get snp info
a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]



snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          snp.id=info[chr!="Scaffold_3", snp.id],
                          sample.id=samps,
              
                          autosome.only=FALSE)

snpset.X<- snpgdsLDpruning(genofile,
                           ld.threshold=0.2,
                           slide.max.bp = 5000,
                           snp.id=info[chr=="Scaffold_3", snp.id],
                           autosome.only=FALSE,
                           sample.id=samps[samps%in%metadata[assigned_sex=="F", sample.id]])

snpset.id <-unlist(snpset)
snpset.id.X <-unlist(snpset.X)


pca.by.chr<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[assigned_sex=="F", sample.id]],
                           num.thread=10)
  }else{
  pca.temp2 <- snpgdsPCA(genofile ,
                         snp.id=snps.to.use, 
                         sample.id=samps,
                         autosome.only=FALSE, 
                         num.thread=10)}
  pca.return <- data.table(pc1=pca.temp2$eigenvect[,1],
                        pc2=pca.temp2$eigenvect[,2],
                        pc3=pca.temp2$eigenvect[,3],
                        pc4=pca.temp2$eigenvect[,4],
                        sample.id=pca.temp2$sample, 
                        chr=i)
  return(pca.return)
}

pca.by.chr<-rbindlist(pca.by.chr)
pca.by.chr<-merge(metadata, pca.by.chr, by="sample.id")


fwrite(pca.by.chr, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_bychr_mac3_noinv.csv")

pca.by.chr[,continent:=factor(continent, levels=c("Africa", "N. America", "S. America", "Hawaii"))]   
pca.by.chr<-pca.by.chr[order(continent)]

pc.chr.all<-ggplot(pca.by.chr, aes(x=pc1, y=pc2, fill=continent))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x="PC1", y="PC2")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 


#####################################
 ##North America by chromosome
#####################################


snpset.na <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. America", sample.id]],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=(3/(2*length(samps[samps%in%metadata[continent=="N. America", sample.id]]))))


snpset.na.X <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. America"&assigned_sex=="F", sample.id]],
                             snp.id=info[chr=="Scaffold_3", snp.id],
                             maf=(3/(2*length(samps[samps%in%metadata[continent=="N. America", sample.id]]))))


snpset.id.na <-unlist(snpset.na)
snpset.id.na.X <-unlist(snpset.na.X)


pca.by.chr.na<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id.na&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.na.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[continent=="N. America"&assigned_sex=="F", sample.id]],
                           num.thread=10)
  }else{
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snps.to.use, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[continent=="N. America", sample.id]],
                           num.thread=10)}
  
  pca.return <- data.table(pc1=pca.temp2$eigenvect[,1],
                           pc2=pca.temp2$eigenvect[,2],
                           pc3=pca.temp2$eigenvect[,3],
                           pc4=pca.temp2$eigenvect[,4],
                           sample.id=pca.temp2$sample, 
                           chr=i)
  return(pca.return)
}

pca.by.chr.na<-rbindlist(pca.by.chr.na)
pca.by.chr.na<-merge(metadata, pca.by.chr.na, by="sample.id")

fwrite(pca.by.chr.na, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_bychr_northamerica_mac3_noinv.csv")

pca.by.chr.na[,loc.spec:=factor(loc.spec, levels=c("Northeast", "VA-CM", "VA-HPO", "TN", "NC", "FL"))]
pca.by.chr.na<-pca.by.chr.na[order(loc.spec)]
  
pc.chr.na<-ggplot(pca.by.chr.na, aes(x=pc1, y=pc2, fill=loc.spec))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x="PC1", y="PC2", fill="Location")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


foreach(i=c(1:5))%do%{
  print(i)
  print("PC1")
  print(summary(aov(pca.by.chr.na[chr==i]$pc1~pca.by.chr.na[chr==i]$loc.spec)))
  print("PC2")
  print(summary(aov(pca.by.chr.na[chr==i]$pc2~pca.by.chr.na[chr==i]$loc.spec)))
}

TukeyHSD(aov(pca.by.chr.na[chr==3]$pc1~pca.by.chr.na[chr==3]$loc.spec))
TukeyHSD(aov(pca.by.chr.na[chr==5]$pc2~pca.by.chr.na[chr==5]$loc.spec))


#####################################
#CM by chromosome
#####################################

snpset.cm <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=3/(2*length(samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]])))


snpset.cm.X <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"&assigned_sex=="F", sample.id]],
                             snp.id=info[chr=="Scaffold_3", snp.id],
                             maf=3/(2*length(samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]])))



snpset.id.cm <-unlist(snpset.cm)
snpset.id.cm.X <-unlist(snpset.cm.X)


pca.by.chr.cm<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id.cm&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.cm.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"&assigned_sex=="F", sample.id]],
                           num.thread=10)
  }else{
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snps.to.use, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]],
                           num.thread=10)}
  
  pca.return <- data.table(pc1=pca.temp2$eigenvect[,1],
                           pc2=pca.temp2$eigenvect[,2],
                           pc3=pca.temp2$eigenvect[,3],
                           pc4=pca.temp2$eigenvect[,4],
                           sample.id=pca.temp2$sample, 
                           chr=i)
  return(pca.return)
}

pca.by.chr.cm<-rbindlist(pca.by.chr.cm)
pca.by.chr.cm<-merge(metadata, pca.by.chr.cm, by="sample.id")

fwrite(pca.by.chr.cm, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_CMonly_mac3_noinv.csv")


pc.chr.cm<-ggplot(pca.by.chr.cm, aes(x=pc1, y=pc2, fill=as.factor(Year)))+
  geom_point(color="white", pch=21, size=4, alpha=0.75)+
  labs(x="PC1", y="PC2", fill="Year")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_fill_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



jpeg("/scratch/perickso/private/ind_seq/Figures/FigureS5_PCAbychr_mac3_noinv.jpg", height=10, width=12, units="in", res=300)
plot_grid(pc.chr.all, pc.chr.na, pc.chr.cm, nrow=3, labels=c("A", "B", "C"), align="v", axis="lr")
dev.off()

foreach(i=c(1:5))%do%{
  print(i)
  print("PC1")
  print(summary(aov(pca.by.chr.cm[chr==i]$pc1~as.factor(pca.by.chr.cm[chr==i]$Year))))
  print("PC2")
  print(summary(aov(pca.by.chr.cm[chr==i]$pc2~as.factor(pca.by.chr.cm[chr==i]$Year))))
  print("PC3")
  print(summary(aov(pca.by.chr.cm[chr==i]$pc3~as.factor(pca.by.chr.cm[chr==i]$Year))))
  print("PC4")
  print(summary(aov(pca.by.chr.cm[chr==i]$pc4~as.factor(pca.by.chr.cm[chr==i]$Year))))
}

TukeyHSD(aov(pca.by.chr.cm[chr==4]$pc1~as.factor(pca.by.chr.cm[chr==4]$Year)))
TukeyHSD(aov(pca.by.chr.cm[chr==4]$pc2~as.factor(pca.by.chr.cm[chr==4]$Year)))