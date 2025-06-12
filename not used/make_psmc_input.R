library(data.table)
library(foreach)
pops<-fread("/scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt", header=F)$V1

input<-foreach(p=pops, .combine="rbind")%do%{
  subpop<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/phlash/psmcfa/", p, ".samps.smc.txt"), header=F)$V1
  return(data.table(pop=p,
                    subpop=subpop))

}

write.table(input, file="/scratch/perickso/private/ind_seq/popgen/phlash/psmc_validation_input.txt", quote=F, sep="\t", row.names=F, col.names=F)