library(gdsfmt)
library(SNPRelate)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(ggsci)
library(ggpubfigs)

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)

females<-metadata[assigned_sex=="F", sample.id]

metadata[loc.spec=="FL", fst.test:="FL"]
metadata[(Location=="VA-CM"|Location=="VA-HPO")&Season!="late", fst.test:="VA"]
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 position=a$pos,
                 freq=a$afreq)

autosomes<-info[chr!="Scaffold_3", snp.id]

X<-info[chr=="Scaffold_3", snp.id]




#calculate pairwise fst

    fst.X<-snpgdsFst(geno, 
                   population=metadata[assigned_sex=="F"&!is.na(fst.test),as.factor(fst.test)],
                   sample.id=metadata[assigned_sex=="F"&!is.na(fst.test),sample.id],
                   autosome.only=F,
                   remove.monosnp=T,
                   maf=.01,
                   missing.rate=0.1,
                   with.id=T, 
                   snp.id=X)
    
    
    fst.auto<-snpgdsFst(geno, 
                     population=metadata[!is.na(fst.test),as.factor(fst.test)],
                     sample.id=metadata[!is.na(fst.test),sample.id],
                     autosome.only=F,
                     remove.monosnp=T,
                     maf=.01,
                     missing.rate=0.1,
                     with.id=T, 
                     snp.id=autosomes)
    z<-data.table(snp.id=c(fst.auto$snp.id,fst.X$snp.id),
                fst.snp=c(fst.auto$FstSNP, fst.X$FstSNP))
    
    
    z[,fst.quant:=frank(fst.snp)/(length(fst.snp)+1)]
z<-merge(z, a, by="snp.id")

save(z, file="/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")

#calculate allele frequency difference (VA vs FL)

FL.auto<-snpgdsSNPRateFreq(geno, 
                           sample.id=metadata[fst.test=="FL",sample.id], 
                           snp.id=autosomes,
                           with.snp.id=T)

FL.X<-snpgdsSNPRateFreq(geno, 
                        sample.id=metadata[fst.test=="FL"&assigned_sex=="F",sample.id], 
                        snp.id=X,
                        with.snp.id=T)


VA.auto<-snpgdsSNPRateFreq(geno, 
                           sample.id=metadata[fst.test=="VA",sample.id], 
                           snp.id=autosomes,
                           with.snp.id=T)

VA.X<-snpgdsSNPRateFreq(geno, 
                        sample.id=metadata[fst.test=="VA"&assigned_sex=="F",sample.id], 
                        snp.id=X,
                        with.snp.id=T)

af_diff=data.table(snp.id=c(FL.auto$snp.id, FL.X$snp.id),
                   freq.FL=c(FL.auto$AlleleFreq, FL.X$AlleleFreq),
                   freq.VA=c(VA.auto$AlleleFreq, VA.X$AlleleFreq))

af_diff[,diff:=abs(freq.FL-freq.VA)]
af_diff<-merge(af_diff, info, by="snp.id")
save(af_diff, file="/scratch/perickso/private/ind_seq/popgen/allelefreqdiff_VAvsFL.Rdata")


#plot for poster:
jpeg("/scratch/perickso/private/ind_seq/popgen/plots/FL-VA_FST_manhattan.jpeg", height=3.5, width=9, res=300, units="in")
ggplot(z)+geom_point(aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_uchicago(palette="dark")+
  labs(x=NULL, y="Florida-Virginia FST")+
  guides(color="none")
dev.off()



#for haplotype plot (in another script)
snps.to.plot<-z[chromosome=="Scaffold_3"&fst.snp>.4&position<1200000]$position

#plot for paper:
jpeg("/scratch/perickso/private/ind_seq/Figures/FL-VA_FST_manhattan.jpeg", height=3.5, width=9, res=300, units="in")
ggplot(z)+geom_point(aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("bright_seven"))+
  labs(x=NULL, y="Florida-Virginia FST")+
  guides(color="none")
dev.off()


#test FST for early vs late in virginia:


metadata[(Location=="VA-CM"|Location=="VA-HPO")&Season!="late", fst.test2:="early"]
metadata[(Location=="VA-CM"|Location=="VA-HPO")&Season=="late", fst.test2:="late"]


fst.X<-snpgdsFst(geno, 
                 population=metadata[assigned_sex=="F"&!is.na(fst.test2),as.factor(fst.test2)],
                 sample.id=metadata[assigned_sex=="F"&!is.na(fst.test2),sample.id],
                 autosome.only=F,
                 remove.monosnp=T,
                 maf=.01,
                 missing.rate=0.1,
                 with.id=T, 
                 snp.id=X)


fst.auto<-snpgdsFst(geno, 
                    population=metadata[!is.na(fst.test2),as.factor(fst.test2)],
                    sample.id=metadata[!is.na(fst.test2),sample.id],
                    autosome.only=F,
                    remove.monosnp=T,
                    maf=.01,
                    missing.rate=0.1,
                    with.id=T, 
                    snp.id=autosomes)
z<-data.table(snp.id=c(fst.auto$snp.id,fst.X$snp.id),
              fst.snp=c(fst.auto$FstSNP, fst.X$FstSNP))


z[,fst.quant:=frank(fst.snp)/(length(fst.snp)+1)]
z<-merge(z, a, by="snp.id")

save(z, file="/scratch/perickso/private/ind_seq/popgen/FST_VAearlyVSlate.Rdata")

#plot for paper: 
jpeg("/scratch/perickso/private/ind_seq/Figures/VAearlyVSlateFST_manhattan.jpeg", height=6, width=8, res=300, units="in")

ggplot(z)+geom_point(aes(x=snp.id, y=fst.snp, color=chromosome))+
  scale_color_manual(values = friendly_pal("bright_seven"))+
  labs(x="SNP ID", y="Virginia early-late FST")+
  guides(color="none")
dev.off()




#loop to look at each population vs Florida

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
metadata[,group:=paste(Location, Year, Season, sep="_")]
groups<-unique(metadata$group)[c(17, 18, 20:26)]

a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos,
                 freq=a$afreq)

autosomes<-info[chr!="Scaffold_3", snp.id]

X<-info[chr=="Scaffold_3", snp.id]

fst.by.pop<-foreach(pop=groups)%do%{
  metadata[,fst.test:=rep(NA, nrow(metadata))]
  metadata[,fst.test:=as.character(fst.test)]
  metadata[loc.spec=="FL", fst.test:="FL"]
  metadata[group==pop, fst.test:=group]
  
  
  #calculate pairwise fst
  
  fst.X<-snpgdsFst(geno, 
                   population=metadata[assigned_sex=="F"&!is.na(fst.test),as.factor(fst.test)],
                   sample.id=metadata[assigned_sex=="F"&!is.na(fst.test),sample.id],
                   autosome.only=F,
                   remove.monosnp=T,
                   maf=.01,
                   missing.rate=0.1,
                   with.id=T, 
                   snp.id=X)
  
  
  fst.auto<-snpgdsFst(geno, 
                      population=metadata[!is.na(fst.test),as.factor(fst.test)],
                      sample.id=metadata[!is.na(fst.test),sample.id],
                      autosome.only=F,
                      remove.monosnp=T,
                      maf=.01,
                      missing.rate=0.1,
                      with.id=T, 
                      snp.id=autosomes)
  z<-data.table(snp.id=c(fst.auto$snp.id,fst.X$snp.id),
                fst.snp=c(fst.auto$FstSNP, fst.X$FstSNP))
  
  
  z[,fst.quant:=frank(fst.snp)/(length(fst.snp)+1)]
  z<-merge(z, a, by="snp.id")
  z[,focal.pop:=pop]
  return(z)
}

fst.by.pop=rbindlist(fst.by.pop)


jpeg("/scratch/perickso/private/ind_seq/Figures/VAvsFL_FST_singlepop.jpeg", height=10, width=10, res=300, units="in")

ggplot(fst.by.pop)+geom_point(aes(x=snp.id, y=fst.snp, color=chromosome))+
  geom_point(data=fst.by.pop[position%in%snps.to.plot&chromosome=="Scaffold_3"], aes(x=snp.id, y=fst.snp), color="#EE6677")+
  scale_color_manual(values = friendly_pal("bright_seven"))+
  labs(x="SNP ID", y="FST")+
  guides(color="none")+
  scale_x_continuous(breaks=c(0, 1000000, 2000000))+
  facet_wrap(~focal.pop)

dev.off()
