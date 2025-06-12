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
metadata[(Location=="VA-CM"|Location=="VA-HPO"), fst.test:="VA"]
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
                 snp.id=X,
                 method = "W&C84")


fst.auto<-snpgdsFst(geno, 
                    population=metadata[!is.na(fst.test),as.factor(fst.test)],
                    sample.id=metadata[!is.na(fst.test),sample.id],
                    autosome.only=F,
                    remove.monosnp=T,
                    maf=.01,
                    missing.rate=0.1,
                    with.id=T, 
                    snp.id=autosomes,
                    method = "W&C84")
z<-data.table(snp.id=c(fst.auto$snp.id,fst.X$snp.id),
              fst.snp=c(fst.auto$FstSNP, fst.X$FstSNP))


z[,fst.quant:=frank(fst.snp)/(length(fst.snp)+1)]
z<-merge(z, a, by="snp.id")

save(z, file="/scratch/perickso/private/ind_seq/popgen/FST_VA_ALLvsFL.Rdata")


###################
#virginia vs africa
###################


geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)

females<-metadata[assigned_sex=="F", sample.id]

metadata[continent=="Africa", fst.test:="Africa"]
metadata[(Location=="VA-CM"|Location=="VA-HPO"), fst.test:="VA"]
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
                 snp.id=X,
                 method = "W&C84")


fst.auto<-snpgdsFst(geno, 
                    population=metadata[!is.na(fst.test),as.factor(fst.test)],
                    sample.id=metadata[!is.na(fst.test),sample.id],
                    autosome.only=F,
                    remove.monosnp=T,
                    maf=.01,
                    missing.rate=0.1,
                    with.id=T, 
                    snp.id=autosomes,
                    method = "W&C84")
z<-data.table(snp.id=c(fst.auto$snp.id,fst.X$snp.id),
              fst.snp=c(fst.auto$FstSNP, fst.X$FstSNP))


z[,fst.quant:=frank(fst.snp)/(length(fst.snp)+1)]
z<-merge(z, a, by="snp.id")

save(z, file="/scratch/perickso/private/ind_seq/popgen/FST_VAvsAfrica.Rdata")


#get allele frequencies for each SNP in each population

FL.auto.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="FL", sample.id], snp.id=autosomes, with.id=T)
FL.X.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="FL"&assigned_sex=="F", sample.id], snp.id=X)
VA.auto.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="VA", sample.id], snp.id=autosomes)
VA.X.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="VA"&assigned_sex=="F", sample.id], snp.id=X)
Af.auto.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="Africa", sample.id], snp.id=autosomes)
Af.X.freq<-snpgdsSNPRateFreq(geno,sample.id=metadata[fst.test=="Africa"&assigned_sex=="F", sample.id], snp.id=X)


freqs<-data.table(snp.id=c(autosomes, X),
                  FL.freq=c(FL.auto.freq$AlleleFreq, FL.X.freq$AlleleFreq),
                  VA.freq=c(VA.auto.freq$AlleleFreq, VA.X.freq$AlleleFreq),
                  Africa.freq=c(Af.auto.freq$AlleleFreq, Af.X.freq$AlleleFreq))

freqs<-merge(freqs, info, by="snp.id")

write.csv(freqs, file="/scratch/perickso/private/ind_seq/popgen/allelefreqsbySNP.csv")
