library(data.table)

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated.csv", header=T)
metadata<-metadata[Location!="INBRED"]
metadata[Location=="CM", Location:="VA-CM"]
metadata[Location=="HPO", Location:="VA-HPO"]
metadata[Location=="", Location:=tstrsplit(tstrsplit(long_name, split="zind-")[[2]], split="_")[[1]]]
metadata[Location=="kenya", Location:="Kenya"]
metadata[Location%in%c("SenegalDesert", "Kenya", "SaoTome", "SenegalForest", "Zambia" ), continent:="Africa"]
metadata[Location=="HI", continent:="Hawaii"]
metadata[Location=="Colombia", continent:="SouthAmerica"]
metadata[is.na(continent), continent:="NorthAmerica"]
setnames(metadata, "sample.name", "sample.id")
metadata[,group:=paste(Location, Year, Season, sep="_")]
metadata[,loc.spec:=Location]
metadata[Location=="TN", loc.spec:="TN"]
metadata[Location=="NC", loc.spec:="NC"]
metadata[Location=="MIA", loc.spec:="FL"]
metadata[Location%in%c("NJ", "PA", "NY"), loc.spec:="Northeast"]
metadata[continent=="Africa", loc.spec:="Africa"]
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
metadata[,group:=paste(Location, Year, Season,  sep="_")]

orig.seq<-metadata[,.(n.seq=.N), .(group, Location, Season,Year, date)]

metadata2<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
act.seq<-metadata2[,.(n.actual=.N), .(group)]

final<-merge(orig.seq, act.seq, by="group")

write.csv(final, "/scratch/perickso/private/ind_seq/sample_sequencing_data.csv", quote=F, row.names=F)