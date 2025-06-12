library(data.table)
library(SNPRelate)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gds", allow.fork=T)
samps.in.geno <- read.gdsn(index.gdsn(geno, "sample.id")) 

females<-metadata[assigned_sex=="F" & sample.id %in% samps.in.geno,.(sample.id)]
females[,fam:=sample.id]


write.table(females, file="/scratch/perickso/private/ind_seq/popgen/admix.female", quote = F, row.names=F, col.names=F)