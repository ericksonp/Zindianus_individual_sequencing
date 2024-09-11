#look at vcf quality

library(data.table)

#first, individual quality

info<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated.csv", header=T)

d<-fread("/scratch/perickso/private/ind_seq/haplotype_calling/updated_Raw_Quality_Statistics.idepth")
names(d)[1]="sample.name"

hist(d$MEAN_DEPTH)

miss<-fread("//scratch/perickso/private/ind_seq/haplotype_calling/updated_Raw_Quality_Statistics.imiss")
names(miss)[1]="sample.name"

hist(miss$F_MISS)

a<-merge(info, merge(d, miss, by="sample.name"), by="sample.name")
#filter individuals missing 10% of genotypes or with mean depth less than 7
a[F_MISS>0.1|MEAN_DEPTH<7]

#write this as a list for the next round of filtering
write.table(a[F_MISS>0.1|MEAN_DEPTH<7, sample.name], "/scratch/perickso/private/ind_seq/haplotype_calling/updated_samples_to_drop.txt", quote=F, row.names=F, sep="\t", col.names=F)

#check for high depth variants
ldepth<-fread("/scratch/perickso/private/ind_seq/haplotype_calling/updated_Raw_Quality_Statistics.ldepth.mean")
hist(ldepth$MEAN_DEPTH)
#remove mean depth > 45 