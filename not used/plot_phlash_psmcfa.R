library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)

pops<-fread("/scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt", header=F)$V1


final<-foreach(pop=pops, .errorhandling="remove")%do%{
  
  p.data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/phlash/", pop, ".autosome.psmcfa.phlash.csv"))
  p.data[,pop:=pop]
  return(p.data)
}
final<-rbindlist(final)
final[,loc:=tstrsplit(pop, split="_")[[1]]]
final[,type:=ifelse(loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"), "native", "introduced")]
final[loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"),type2:="native"]
final[pop=="MIA_2019_June", type2:="introduced-Florida"]
final[grepl("early", pop),type2:="introduced-Virginia"]
final[grepl("mid", pop),type2:="introduced-Virginia"]


final[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
final<-final[order(type2)]
final[,pop:=as.factor(pop)]
final[,pop:=factor(pop, levels=rev(levels(final$pop)))]

phlash.plot<-ggplot()+
  # geom_line(data=folds[type2%in%c("native", "introduced-Florida", "introduced-early")], aes(x=x, y=y, group =interaction(pop,fold), color=as.factor(type2)), size=0.1)+
  geom_line(data=final[type2%in%c("native", "introduced-Florida", "introduced-Virginia")&time<1000], aes(x=time, y=upper97.5, color=type2, group=pop), size=1.5)+  
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  #scale_y_log10(
  #  breaks = scales::trans_breaks("log10", function(x) 10^x),
  #  labels = scales::trans_format("log10", scales::math_format(10^.x))
  #) +
 # annotation_logticks() +
  labs(x="years before present", y="population size", color=NULL)+
  scale_color_manual(values = friendly_pal("contrast_three"))+
  theme(legend.position = c(0.5, 0.8))