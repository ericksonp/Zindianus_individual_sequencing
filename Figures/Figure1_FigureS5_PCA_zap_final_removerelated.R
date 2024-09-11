

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
 #snpgdsVCF2GDS(vcf.fn = "/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz",  
#              out.fn="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds", 
 #             method="biallelic.only")


#load genofile
genofile <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
samps <- read.gdsn(index.gdsn(genofile, "sample.id")) 

#get metadata
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
metadata[continent=="NorthAmerica", continent:="N. Am"]
metadata[continent=="SouthAmerica", continent:="S. Am"]

#get snp info
a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]


#read in unrelated individuals
unrelated<-fread("/scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id")
names(unrelated)="id"

#get informative SNPs excluding X chromosome and filtering at MAF <0.05
snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          sample.id=unrelated$id,
                          snp.id=info[chr!="Scaffold_3", snp.id],
                          autosome.only=FALSE,
                          maf=.05)


#get informative SNPs excluding including X and only females
snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          sample.id=samps[samps%in%unrelated$id&samps%in%metadata[assigned_sex=="F", sample.id]],
                          autosome.only=FALSE,
                          maf=.05)


#turn SNPsets into vectors
snpset.id <-unlist(snpset)


 #calculate PCs

     pca.temp <- snpgdsPCA(genofile 
                           ,snp.id=snpset.id, 
                           autosome.only=FALSE, 
                           num.thread=20, 
                           sample.id=unrelated$id)

     #snpgdsPCASNPLoading(pca.temp, genofile)

     pca.dt <- data.table(pc1=pca.temp$eigenvect[,1],
                          pc2=pca.temp$eigenvect[,2],
                          pc3=pca.temp$eigenvect[,3],
                          pc4=pca.temp$eigenvect[,4],
                          sample.id=pca.temp$sample)
    


 pca.dt=merge(pca.dt, metadata, by="sample.id")
fwrite(pca.dt, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated.csv")
 #plot overall PCA
 pca.percent <- data.table(pc1=pca.temp$eigenval[1]/sum(pca.temp$eigenval, na.rm=T)*100,
                          pc2=pca.temp$eigenval[2]/sum(pca.temp$eigenval, na.rm=T)*100)
                   
         
                  
pa<-ggplot(pca.dt, aes(x=pc1, y=pc2, color=continent))+
  geom_point()+
  labs(x=paste("PC1,", round(pca.percent[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent[1,2], 2), "%", sep=" "))+
 scale_color_manual(values = friendly_pal("ito_seven"))

#calculate stats for North America-Africa PC1 divergence
 t.test(pca.dt[continent%in%c("N. Am", "Africa")]$pc1~pca.dt[continent%in%c("N. Am", "Africa")]$continent)
                

################################
#PCA only North American samples
###############################

#north american samples, autosomes only
#recalcuate ld pruning
snpset.na <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. Am", sample.id] & samps%in%unrelated$id],
                             snp.id=info[chr!="Scaffold_3",snp.id],
                             maf=.05)


snpset.id.na <-unlist(snpset.na)

pca.temp2 <- snpgdsPCA(genofile ,
                       snp.id=snpset.id.na, 
                       autosome.only=FALSE, 
                       sample.id=samps[samps%in%metadata[continent=="N. Am", sample.id] & samps%in%unrelated$id],
                       num.thread=10)

pca.dt2 <- data.table(pc1=pca.temp2$eigenvect[,1],
                      pc2=pca.temp2$eigenvect[,2],
                      pc3=pca.temp2$eigenvect[,3],
                      pc4=pca.temp2$eigenvect[,4],
                      sample.id=pca.temp2$sample)
pca.percent2 <- data.table(pc1=pca.temp2$eigenval[1]/sum(pca.temp2$eigenval, na.rm=T)*100,
                           pc2=pca.temp2$eigenval[2]/sum(pca.temp2$eigenval, na.rm=T)*100,
                           pc3=pca.temp2$eigenval[3]/sum(pca.temp2$eigenval, na.rm=T)*100)


pca.dt2=merge(pca.dt2, metadata, by="sample.id")

fwrite(pca.dt2, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_northamerica.csv")

pb<-ggplot(pca.dt2, aes(x=pc1, y=pc2, color=loc.spec))+geom_point()+
  labs(x=paste("PC1,", round(pca.percent2[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent2[1,2], 2), "%", sep=" "),
       color="Location")+
  scale_color_manual(values = friendly_pal("ito_seven"))

summary(aov(pca.dt2$pc1~pca.dt2$loc.spec))
summary(aov(pca.dt2$pc2~pca.dt2$loc.spec))


########################################
#Carter Mountain by year
########################################


snpset.cm <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]& samps%in%unrelated$id],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=.05)


snpset.id.cm <-unlist(snpset.cm)

pca.temp3 <- snpgdsPCA(genofile ,
                       snp.id=snpset.id.cm, 
                       autosome.only=FALSE, 
                       sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]& samps%in%unrelated$id],
                       num.thread=10)

pca.dt3 <- data.table(pc1=pca.temp3$eigenvect[,1],
                      pc2=pca.temp3$eigenvect[,2],
                      pc3=pca.temp3$eigenvect[,3],
                      pc4=pca.temp3$eigenvect[,4],
                      sample.id=pca.temp3$sample)

pca.percent3 <- data.table(pc1=pca.temp3$eigenval[1]/sum(pca.temp3$eigenval, na.rm=T)*100,
                           pc2=pca.temp3$eigenval[2]/sum(pca.temp3$eigenval, na.rm=T)*100,
                           pc3=pca.temp3$eigenval[3]/sum(pca.temp3$eigenval, na.rm=T)*100)

metadata[Season=="mid", Season:="early"]
pca.dt3=merge(pca.dt3, metadata, by="sample.id")


fwrite(pca.dt3, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_CMonly.csv")

pc<-ggplot(pca.dt3, aes(x=pc1, y=pc2, color=as.factor(Year)))+geom_point()+
  labs(x=paste("PC1,", round(pca.percent3[1,1], 2), "%", sep=" "),
       y=paste("PC2,", round(pca.percent3[1,2], 2), "%", sep=" "),
       color="Year")+
  scale_color_manual(values = friendly_pal("ito_seven"))

#stats
summary(aov(pca.dt3$pc1~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc2~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc3~as.factor(pca.dt3$Year)))
summary(aov(pca.dt3$pc4~as.factor(pca.dt3$Year)))


pdf("/scratch/perickso/private/ind_seq/Figures/Figure_1_PCA_unrelated.pdf", height=4, width=12)
plot_grid(pa, pb, pc, nrow=1, labels=c("A", "B", "C"))
dev.off()





##########################
### FOR SUPPLEMENT: repeat all analyses by chromosome
#do PCA by individual chromosomes; only use females for scaffold 3
#########################



snpset <- snpgdsLDpruning(genofile,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          snp.id=info[chr!="Scaffold_3", snp.id],
                          sample.id=unrelated$id,
              
                          autosome.only=FALSE,
                          maf=.05)

snpset.X<- snpgdsLDpruning(genofile,
                           ld.threshold=0.2,
                           slide.max.bp = 5000,
                           snp.id=info[chr=="Scaffold_3", snp.id],
                           autosome.only=FALSE,
                           sample.id=samps[samps%in%metadata[assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                           maf=.05)

snpset.id <-unlist(snpset)
snpset.id.X <-unlist(snpset.X)


pca.by.chr<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                           num.thread=10)
  }else{
  pca.temp2 <- snpgdsPCA(genofile ,
                         snp.id=snps.to.use, 
                         sample.id=samps[samps%in%unrelated$id],
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


fwrite(pca.by.chr, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr.csv")



pc.chr.all<-ggplot(pca.by.chr, aes(x=pc1, y=pc2, color=continent))+
  geom_point()+labs(x="PC1", y="PC2")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 


pc34.chr.all<-ggplot(pca.by.chr, aes(x=pc3, y=pc4, color=continent))+
  geom_point()+labs(x="PC3", y="PC4")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#####################################
 ##North America by chromosome
#####################################


snpset.na <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. Am", sample.id]& samps%in%unrelated$id],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=.05)

snpset.na.X <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[continent=="N. Am"&assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                             snp.id=info[chr=="Scaffold_3", snp.id],
                             maf=.05)

snpset.id.na <-unlist(snpset.na)
snpset.id.na.X <-unlist(snpset.na.X)


pca.by.chr.na<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id.na&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.na.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[continent=="N. Am"&assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                           num.thread=10)
  }else{
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snps.to.use, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[continent=="N. Am", sample.id]& samps%in%unrelated$id],
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

fwrite(pca.by.chr.na, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_northamerica.csv")
foreach(i=c(1:5))%do%{
  print(i)
  print("PC1")
  print(summary(aov(pca.by.chr.na[chr==i]$pc1~pca.by.chr.na[chr==i]$loc.spec)))
  print("PC2")
  print(summary(aov(pca.by.chr.na[chr==i]$pc2~pca.by.chr.na[chr==i]$loc.spec)))
}

TukeyHSD(aov(pca.by.chr.na[chr==3]$pc1~pca.by.chr.na[chr==3]$loc.spec))
TukeyHSD(aov(pca.by.chr.na[chr==5]$pc2~pca.by.chr.na[chr==5]$loc.spec))

  
pc.chr.na<-ggplot(pca.by.chr.na, aes(x=pc1, y=pc2, color=loc.spec))+
  geom_point()+
  labs(x="PC1", y="PC2", color="Location")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pc34.chr.na<-ggplot(pca.by.chr.na, aes(x=pc3, y=pc4, color=loc.spec))+
  geom_point()+
  labs(x="PC3", y="PC4", color="Location")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




#####################################
#CM by chromosome
#####################################

snpset.cm <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]& samps%in%unrelated$id],
                             snp.id=info[chr!="Scaffold_3", snp.id],
                             maf=.05)

snpset.cm.X <- snpgdsLDpruning(genofile,
                             ld.threshold=0.2,
                             slide.max.bp = 5000,
                             autosome.only=FALSE,
                             sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"&assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                             snp.id=info[chr=="Scaffold_3", snp.id],
                             maf=.05)


snpset.id.cm <-unlist(snpset.cm)
snpset.id.cm.X <-unlist(snpset.cm.X)


pca.by.chr.cm<-foreach(i=c(1:5))%do%{
  print("chr ", i)
  snps.to.use<-info[snp.id%in%snpset.id.cm&chr==paste0("Scaffold_", i), snp.id]
  if(i==3){
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snpset.id.cm.X, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[loc.spec=="VA-CM"&assigned_sex=="F", sample.id]& samps%in%unrelated$id],
                           num.thread=10)
  }else{
    pca.temp2 <- snpgdsPCA(genofile ,
                           snp.id=snps.to.use, 
                           autosome.only=FALSE, 
                           sample.id=samps[samps%in%metadata[loc.spec=="VA-CM", sample.id]& samps%in%unrelated$id],
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

fwrite(pca.by.chr.cm, file= "/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_CMonly.csv")

pca.by.chr.cm<-fread("/scratch/perickso/private/ind_seq/popgen/PCA_all_unrelated_bychr_CMonly.csv")

pc.chr.cm<-ggplot(pca.by.chr.cm, aes(x=pc1, y=pc2, color=as.factor(Year)))+
  geom_point()+labs(x="PC1", y="PC2", color="Year")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pc34.chr.cm<-ggplot(pca.by.chr.cm, aes(x=pc3, y=pc4, color=as.factor(Year)))+
  geom_point()+labs(x="PC3", y="PC4", color="Year")+
  facet_wrap(~chr, scales="free", ncol=5)+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

jpeg("/scratch/perickso/private/ind_seq/Figures/FigureS5_PCAbychr_unrelated.jpg", height=10, width=12, units="in", res=300)
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