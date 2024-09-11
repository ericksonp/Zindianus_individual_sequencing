library(data.table)
unrelated<-fread("/scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id")
names(unrelated)="id"

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)

unrelated<-unrelated[id%in%metadata[assigned_sex=="F", sample.id]]
unrelated[,fam:=id]


write.table(unrelated, file="/scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id.female", quote = F, row.names=F, col.names=F)