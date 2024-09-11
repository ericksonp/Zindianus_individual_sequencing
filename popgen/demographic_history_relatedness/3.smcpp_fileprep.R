library(data.table)
library(foreach)

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv")
write.table(unique(metadata$group), file="/scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt", quote=F, sep="\t", row.names=F, col.names=F)

foreach(pop=unique(metadata$group))%do%{
  print(pop)
  samps<-metadata[group==pop, sample.id]
  write.table(samps, file=paste0("/scratch/perickso/private/ind_seq/popgen/smcpp/", pop, ".samps.smc.txt"), quote=F, sep="\t", row.names=F, col.names=F)
}


#make other populations not just based on collection location and year:

#NC

write.table(metadata[Location=="NC", sample.id], file="/scratch/perickso/private/ind_seq/popgen/smcpp/NC.samps.smc.txt", quote=F, sep="\t", row.names=F, col.names=F)
write.table(metadata[loc.spec=="Northeast", sample.id], file="/scratch/perickso/private/ind_seq/popgen/smcpp/Northeast.samps.smc.txt", quote=F, sep="\t", row.names=F, col.names=F)
write.table(metadata[loc.spec=="FL", sample.id], file="/scratch/perickso/private/ind_seq/popgen/smcpp/FL.samps.smc.txt", quote=F, sep="\t", row.names=F, col.names=F)

#add 3 samps to pops_for_smcpp.txt with nano and run slurm with jobs 27-29