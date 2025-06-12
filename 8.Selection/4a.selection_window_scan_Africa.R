library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)
library(gdsfmt)
library(SNPRelate)
library(lattice)
library(tidyr)
library(stringr)
library(lubridate)
library(viridis)
library(ggpubfigs)
library(rehh)
library(scales)
library(doMC)
registerDoMC(50)


#FST (loads in as z with SNPids)
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsAfrica.Rdata")
z[,old.snp.id:=snp.id]
fst.99<-quantile(z$fst.snp, 0.99, na.rm=T)
z[,q99:=fst.snp>=fst.99]

#IHS

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VAall<-as.data.table(wgscan.ihs$ihs)

geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)

a<-snpgdsSNPList(geno)
info<-data.table(old.snp.id=a$snp.id,
                 CHR=a$chromosome,
                 POSITION=a$pos,
                 freq=a$afreq)
ihs.VAall<-merge(ihs.VAall, info, by=c("CHR", "POSITION"))
ihs.99<-quantile(ihs.VAall$IHS, 0.99, na.rm=T)
ihs.VAall[,q99:=IHS>=ihs.99]






#Bayescan

bp.va.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/updated_data_for_manhattan.txt")
bp.va.af[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]

setnames(info, "POSITION", "pos")
info[,Scaffold:=as.integer(tstrsplit(CHR, split="_")[[2]])]
bp.va.af<-merge(bp.va.af, info, by=c("pos", "Scaffold"))
bp.99<-quantile(bp.va.af$M_XtX, 0.99, na.rm=T)
bp.va.af[,q99:=M_XtX>=bp.99]

#try looking for individual SNPs by merging all together-have to get same column names here

### get everything the same

setnames(z, c("chromosome", "position" ), c("chr", "pos"))
setnames(z, "q99", "fstq99")
setnames(ihs.VAall, c("CHR", "POSITION"), c("chr", "pos"))
setnames(ihs.VAall, "q99", "ihsq99" )
setnames(bp.va.af, "CHR", "chr")
setnames(bp.va.af, "q99", "bpq99")

#merge
alldata<-merge(z, ihs.VAall, by=c("chr", "pos" , "old.snp.id"), all=T)
alldata<-merge(alldata, bp.va.af, by=c("chr", "pos", "old.snp.id"), all=T)

#find snps shared across two or three tests
alldata[,all.3:=(fstq99==T & ihsq99==T & bpq99==T)]
alldata[,bp.fst:=fstq99==T&bpq99==T]
alldata[,bp.ihs:=bpq99==T&ihsq99==T]
alldata[,ihs.fst:=ihsq99==T&fstq99==T]


#look at windows 
win.bp <- 1000
step.bp <- 500

wins<-foreach(chr.i=c(1:5),
              .combine="rbind", 
              .errorhandling="remove")%dopar%{
                
                tmp <- info[CHR==paste0("Scaffold_",chr.i)]
                data.table(CHR=chr.i,
                           start=seq(from=1, to=max(tmp$pos)-(win.bp), by=step.bp),
                           end=seq(from=0, to=max(tmp$pos)-(win.bp), by=step.bp) + win.bp)
              }

wins[,index:=c(1:nrow(wins))]
setkey(wins,index)

setkey(z, chr, pos)
setkey(ihs.VAall, chr, pos)
setkey(bp.va.af, chr, pos)
window.sum<-foreach(window=wins$index, .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  if(window%%10000==0){
                    print(window)
                  }
                  fst.tmp<-z[chr==paste0("Scaffold_", wins[J(window),CHR]) & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  fst.count<-sum(fst.tmp$fstq99, na.rm=T)
                  ihs.tmp<-ihs.VAall[chr==paste0("Scaffold_", wins[J(window),CHR]) & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  ihs.count<-sum(ihs.tmp$ihsq99, na.rm=T)
                  bp.tmp<-bp.va.af[chr==paste0("Scaffold_", wins[J(window),CHR]) & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  bp.count<-sum(bp.tmp$bpq99, na.rm=T)
            
                  data.table(index=window,
                             fst.count=fst.count,
                             n.fst=nrow(fst.tmp),
                             ihs.count=ihs.count,
                             n.ihs=nrow(ihs.tmp),
                             bp.count=bp.count,
                             n.bp=nrow(bp.tmp))
                            
                             
                }

window.sum<-merge(window.sum, wins, by="index")
window.sum[,all.3:=fst.count>0&ihs.count>0&bp.count>0]
window.sum[,bp.fst:=fst.count>0&bp.count>0]
window.sum[,bp.ihs:=bp.count>0&ihs.count>0]
window.sum[,ihs.fst:=ihs.count>0&fst.count>0]
window.sum[]

#write.csv(window.sum, file="/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_10kbwindows.csv")
write.csv(window.sum, file="/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_Africa.csv")


#cycle through windows to grab SNPs in those windows
setkey(window.sum, index)
snptags<-foreach (window=window.sum[all.3==T,index], .combine="rbind", .errorhandling="remove")%dopar%{
  tmp.chr=window.sum[window,CHR]
  tmp.start=window.sum[window, start]
  tmp.end=window.sum[window, end]
  snps<-info[Scaffold==tmp.chr&pos>=tmp.start&pos<=tmp.end]
  return(data.table(index=window,
                    chr=paste0("Scaffold_", tmp.chr),
                    first.snp=min(snps$old.snp.id),
                    last.snp=max(snps$old.snp.id)))
  }

write.csv(snptags, file="/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows_Africa_snpstoplot.csv")

