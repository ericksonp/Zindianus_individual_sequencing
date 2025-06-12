library(data.table)
f<-data.table(file_name=list.files("/scratch/perickso/private/ind_seq/haplotype_calling", pattern=".vcf.gz"))
f<-f[!grepl("tbi", file_name)]
f[,sample_name:=tstrsplit(file_name, split="[.]")[[1]]]
#remove blanks
f<-f[sample_name!="ZAP_1_H12"]
f<-f[sample_name!="ZAP_2_H11"]
write.table(f[,.(sample_name, file_name)], "/scratch/perickso/private/ind_seq/Samples_to_haplotype.txt" , quote=F, row.names=F, col.names=F, sep="\t")


#for published list

d<-fread("/scratch/perickso/private/raw_data/published/SRR_Acc_List.txt", header=F)

d[,V2:=paste0(V1, ".raw.g.vcf.gz")]

write.table(d, "/scratch/perickso/private/ind_seq/Samples_to_haplotype_published.txt" , quote=F, row.names=F, col.names=F, sep="\t")
