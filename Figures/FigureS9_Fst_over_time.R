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


#read in unrelated individuals
#unrelated<-fread("/scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id")
#names(unrelated)="id"

#metadata<-metadata[sample.id%in%unrelated$id]
metadata<-metadata[loc.spec=="VA-CM"|loc.spec=="VA-HPO"]

sampling<-unique(metadata[, .(group,date)])
sampling[,date:=mdy(date)]
sampling.avg<-sampling[,.(mean.date=mean(date)), .(group)]
sampling.avg[,pop1:=group]
sampling.avg[,pop2:=group]

#get snp info
a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)
info[,maf:=pmin(freq, 1-freq)]


pop.pairs<-as.data.table(t(combn( unique(metadata$group), 2)))

#calculate pairwise fst
y<-foreach(i=c(1:nrow(pop.pairs)), .errorhandling="remove")%do%{
  pop1=pop.pairs[i, V1]
  pop2=pop.pairs[i,V2]
  submeta<-metadata[group==pop1|group==pop2]
  fst<-snpgdsFst(genofile, 
                 population=submeta[,as.factor(group)],
                 sample.id=submeta[,sample.id],
                 autosome.only=F,
                 remove.monosnp=T,
                 snp.id=info[chr!="Scaffold_3", snp.id],
                 maf=.05,
                 method="W&C84")
  return(data.table(pop1=pop1,
                    pop2=pop2,
                    Fst=fst$Fst,
                    MeanFst=fst$MeanFst))
}

y<-rbindlist(y)

y<-merge(y, sampling.avg[,.(pop1,mean.date)], by="pop1")
setnames(y, "mean.date", "pop1.date")

y<-merge(y, sampling.avg[,.(pop2,mean.date)], by="pop2")
setnames(y, "mean.date", "pop2.date")

y[,loc1:=tstrsplit(pop1, split="_")[[1]]]
y[,loc2:=tstrsplit(pop2, split="_")[[1]]]

y[,timediff:=abs(pop1.date-pop2.date)]

ggplot(y[loc1=="VA-CM"&loc2=="VA-CM"])+geom_point(aes(x=timediff, y=Fst))+
  geom_smooth(aes(x=timediff, y=Fst), method="lm")+
  labs(x="Days between collections",
       y=expression("F"[ST]))

summary(lm(Fst~timediff, data=y[loc1=="VA-CM"&loc2=="VA-CM"])) #P=0.32

y[,weird:=ifelse(pop1=="VA-CM_2018_early"|pop2=="VA-CM_2018_early", "yes", "no")]

jpeg("/scratch/perickso/private/ind_seq/Figures/Fst_over_time.jpg", height=5, width=5, units="in", res=300)
ggplot(y[loc1=="VA-CM"&loc2=="VA-CM"])+geom_point(aes(x=timediff, y=Fst, color=as.factor(weird)))+
  geom_smooth(aes(x=timediff, y=Fst, color=as.factor(weird)), method="lm", se=F)+
  labs(x="Days between collections",
       y=expression("F"[ST]),
       color="comparing to\n2018 early")+
  scale_color_manual(values = friendly_pal("contrast_three"))
dev.off()
  
summary(lm(Fst~timediff*weird, data=y[loc1=="VA-CM"&loc2=="VA-CM"])) #P=0.32
