library(data.table)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
treemix<-metadata[Location%in%c("VA-CM", "VA-HPO", "MIA")]
cluster<-data.table(id=treemix$sample.id,
                    id2=treemix$sample.id,
                    group=treemix$group)

write.table(cluster[,id], file="/scratch/perickso/private/ind_seq/popgen/treemix/treemix_samps.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(cluster, file="/scratch/perickso/private/ind_seq/popgen/treemix/treemix_cluster.txt", quote=F, col.names=F, row.names=F, sep="\t")


#try just CM populations and MIA, group by year rather than season

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
treemix<-metadata[Location%in%c("VA-CM", "MIA")]
cluster<-data.table(id=treemix$sample.id,
                    id2=treemix$sample.id,
                    group=paste(treemix$Location, treemix$Year, sep="_"))

write.table(cluster[,id], file="/scratch/perickso/private/ind_seq/popgen/treemix/treemix_samps_CM_MIA.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(cluster, file="/scratch/perickso/private/ind_seq/popgen/treemix/treemix_cluster_CM_MIA.txt", quote=F, col.names=F, row.names=F, sep="\t")