#!/usr/bin/env Rscript

#LD analysis
library(data.table)
library(SNPRelate)
library(foreach)
library(doMC)
registerDoMC(50)


geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gds", allow.fork=T)
samps.in.geno <- read.gdsn(index.gdsn(geno, "sample.id")) 

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
samps<-metadata[continent=="NorthAmerica"&sample.id%in%samps.in.geno, sample.id]

a<-snpgdsSNPList(geno, sample.id=samps)

info<-data.table(chr=a$chromosome,
                 pos=a$position,
                 snp.id=a$snp.id,
                 af=a$afreq)

info[,missing:=snpgdsSNPRateFreq(geno, sample.id=samps)$MissingRate]
info<-info[missing==0&af>0&af<1]

#look at long distance LD decay in random SNPs from swarm

m<-foreach(focal.snp=info$snp.id, .errorhandling="remove")%dopar%{
  if(focal.snp%%10000==0){
    print(focal.snp)
  }
  setkey(info, snp.id)
  focal.chr<-info[J(focal.snp), chr]
  focal.pos<-info[J(focal.snp), pos]
  y<-foreach(dist=c(50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000), .errorhandling = "remove")%do%{
    range<-info[chr==focal.chr&pos>(focal.pos+.95*dist)&pos<(focal.pos+1.05*dist)]
    #print(range)
    focal.snp2<-sample(range$snp.id, 1)
    ld<-snpgdsLDMat(geno, slide=-1,verbose=F,snp.id = c(focal.snp, focal.snp2), method='composite')$LD[2,1]
    
    return(data.table(dist=dist,
                      ld=ld,
                      focal.snp=focal.snp,
                      focal.chr=focal.chr,
                      focal.snp2=focal.snp2,
                      focal.pos=focal.pos,
                      focal.pos2=range[snp.id==focal.snp2, pos]))
  }
  return(rbindlist(y))
  
}
m<-rbindlist(m)

write.table(m, "/scratch/perickso/private/ind_seq/popgen/ld_decay_NorthAmerica_50kb.txt", quote=F, sep="\t", row.names = F)

snpgdsClose(geno)
