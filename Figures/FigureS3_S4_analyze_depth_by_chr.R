library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggsci)
library(ggpubfigs)

z<-foreach(chr=c(1:5))%do%{
  a<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/final_quality_by_chr.", chr, ".idepth"),header=T )
  a[,chrom:=chr]
  return(a)     
  }
z<-rbindlist(z)

z.auto<-z[chr!=3,.(mean.auto.depth=mean(MEAN_DEPTH)), .(INDV)]

z.all<-merge(z.auto, z[chrom==3], by="INDV")

z.all[,sex.to.auto:=MEAN_DEPTH/mean.auto.depth]
setnames(z.all, "INDV", "sample.id")

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", header=T, drop=1)

z.all<-merge(metadata[,.(sample.id, Sex)], z.all, by="sample.id")
z.all[Sex=="", Sex:="Unknown"]
z.all[Sex=="F", Sex:="Female"]
z.all[Sex=="M", Sex:="Male"]


jpeg("/scratch/perickso/private/ind_seq/Figures/sex_to_autosome.jpeg", height=5, width=4, units="in", res=300)
ggplot(z.all)+geom_histogram(aes(x=sex.to.auto))+
  facet_grid(Sex~.)+
  labs(x="Ratio of average sex:autosome SNP depth")+
  geom_vline(xintercept=0.8, linetype="dashed")
dev.off()

z.all[, assigned_sex:=ifelse(sex.to.auto>0.8, "F", "M")]
z.all[Sex!=assigned_sex]


#now plot depth based on population

d<-fread("/scratch/perickso/private/ind_seq/haplotype_calling/updated_Second_Quality_Statistics.idepth", header=T)

setnames(d, "INDV", "sample.id")
d<-merge(metadata,d, by="sample.id")

d[,label:=Location]
d[Season!="", label:=paste(loc.spec, Year, Season, sep="-")]
d[Location=="FL", label:="FL-2016"]
d[Location=="MIA", label:="FL-2019"]

d[,label:=factor(label, levels=c("SaoTome", "SenegalForest", "SenegalDesert", "Kenya", "Zambia", "Colombia", "HI", "FL-2016", "FL-2019", "NC", "TN", "VA-CM-2017-early", "VA-CM-2017-late", "VA-CM-2018-early", "VA-CM-2018-late", "VA-CM-2019-mid", "VA-CM-2019-late", "VA-CM-2020-early", "VA-HPO-2019-mid", "VA-HPO-2020-early", "PA", "NJ", "NY"))]

d[,study:=ifelse(library=="comeault2020"|library=="comeault2021", "Comeault 2020/2021", "this study")]

#plot for paper
jpeg("/scratch/perickso/private/ind_seq/Figures/depth_by_location.jpg", width=8, height=6, units="in", res=300)
ggplot(d)+
  geom_boxplot(aes(x=label, y=MEAN_DEPTH, color=study))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x=NULL, y="Mean SNP depth")+
  scale_color_manual(values = friendly_pal("bright_seven"))
dev.off()

#write.table(z.all[,.(sample.name, sex.to.auto, assigned_sex)], file="/scratch/perickso/private/ind_seq/assigned_sex_by_coverage.txt", quote=F, sep="\t", row.names=F)

