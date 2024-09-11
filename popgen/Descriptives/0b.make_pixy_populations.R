library(data.table)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)

metadata[Location=="VA-CM"|Location=="VA-HPO", pixy:="VA"]
metadata[Location=="MIA"|Location=="FL", pixy:="FL"]
metadata[continent=="Africa", pixy:="Africa"]
metadata[is.na(pixy), pixy:="other"]

write.table(metadata[,.(sample.id, pixy)], "/scratch/perickso/private/ind_seq/popgen/populations_for_pixy.txt", row.names=F, col.names=F, quote=F, sep="\t")

write.table(metadata[assigned_sex=="F",.(sample.id, pixy)], "/scratch/perickso/private/ind_seq/popgen/populations_for_pixy_females.txt", row.names=F, col.names=F, quote=F, sep="\t")

#make list of females only to remove from vcf for analysis of X chromosome

write.table(metadata[assigned_sex=="F", sample.id], "/scratch/perickso/private/ind_seq/popgen/females_for_pixy.txt", row.names=F, col.names=F, quote=F, sep="\t")

seas<-metadata[group%in%unique(metadata$group)[17:26]]

write.table(seas[,.(sample.id, group)], "/scratch/perickso/private/ind_seq/popgen/subpopulations_for_pixy.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(seas[assigned_sex=="F",.(sample.id, group)], "/scratch/perickso/private/ind_seq/popgen/subpopulations_for_pixy_females.txt", row.names=F, col.names=F, quote=F, sep="\t")
