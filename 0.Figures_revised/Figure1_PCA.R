

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


###Make plots #########

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
plot_grid(pa, pb, pc, pap,pbp,pcp, nrow=2, labels=c("a", "b", "c", "a'", "b'", "c'"), align="vh", axis="tblr")
dev.off()


