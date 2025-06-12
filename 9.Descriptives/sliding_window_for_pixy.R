

library(data.table)
library(vcfR)

i=1
read.vcfR(paste0("/scratch/perickso/private/ind_seq/haplotype_calling/Scaffold_",{i},".updated.genotypedSNPs.allsites.raw.recode.vcf.gz"))

v<-as.data.table(getFIX(vcf))
v[,snp_index:=c(1:nrow(v))]

window_size=100
step_size=50


windows<-data.table(CHR=i,
           start=seq(from=1, to=nrow(v)-window_size, by=step_size),
           end=seq(from=1,  to=nrow(v)-window_size, by=step_size) + step_size)


setkey(v, snp_index)

windows[,bed_start:=v[J(start), POS-1]]
windows[,bed_end:=v[J(end), POS-1]]

write.table(windows[,.(CHR, bed_start, bed_end)], paste0("/scratch/perickso/private/ind_seq/popgen/pixy/Scaffold_",i,"_pixy_windows.bed", quote=F, row.names=F, col.names=F, sep="\t")