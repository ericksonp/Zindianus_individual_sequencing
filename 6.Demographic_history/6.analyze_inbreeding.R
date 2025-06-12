library(data.table)
library(ggplot2)
library(beeswarm)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1, header=T)
metadata[,loc.spec:=factor(loc.spec, levels=c("Africa", "Colombia", "FL", "MIA", "NC", "TN", "VA-HPO", "VA-CM", "Northeast", "HI"))]
het<-fread("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.het")
setnames(het, "INDV", "sample.id")

het<-merge(metadata, het, by="sample.id")
#ggplot(het[(loc.spec=="FL"|loc.spec=="VA-CM")&group!="FL_2016_"])+geom_point(aes(x=group, y=F))
summary(aov(F~group, data=het[loc.spec=="VA-CM"]))
