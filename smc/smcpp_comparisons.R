library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(RColorBrewer)
library(ggsci)
library(ggpubfigs)
library(gdsfmt)
library(SNPRelate)
library(lattice)
library(tidyr)
library(stringr)
library(lubridate)
library(viridis)
library(ggpubfigs)

#compare SMCPP runs


pops<-fread("/scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt", header=F)$V1

final<-foreach(pop=pops, .errorhandling="remove")%do%{
  
  p.data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/smcpp/", pop, "/", pop, "_cv_10folds_4knots_10000gen_scaf4_cubic/", pop, "_CV10fold_final_plot.csv"))
  p.data[,pop:=pop]
  return(p.data)
}
final<-rbindlist(final)
final[,loc:=tstrsplit(pop, split="_")[[1]]]
final[,ky:=x/1000]
final[,type:=ifelse(loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"), "native", "introduced")]
final[loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"),type2:="native"]
final[label=="MIA_2019_June", type2:="introduced-Florida"]
final[grepl("early", label),type2:="introduced-Virginia"]
final[grepl("mid", label),type2:="introduced-Virginia"]

final<-final[!grepl("late", label)]
final[is.na(type2), type2:="introduced"]
final[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
final<-final[order(type2)]
final[,label:=as.factor(label)]
final[,label:=factor(label, levels=rev(levels(final$label)))]

folds<-foreach(pop=pops, .errorhandling="remove")%do%{
  
  folds2<-foreach(i=c(0:30), .errorhandling="remove")%do%{
    p.data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/smcpp/", pop, "/", pop, "_cv_10folds_4knots_10000gen_scaf4_cubic/fold", i, "/", pop, "_fold", i, "_final_plot.csv"))
    p.data[,fold:=i]
    p.data[,pop:=pop]
    return(p.data)
    
  }
  return(rbindlist(folds2))
}

folds<-rbindlist(folds)
folds[,loc:=tstrsplit(pop, split="_")[[1]]]
folds[,type:=ifelse(loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"), "native", "introduced")]
folds[loc%in%c("Kenya", "Zambia", "SenegalForest", "SenegalDesert", "SaoTome"),type2:="native"]
folds[label=="MIA_2019_June", type2:="introduced-Florida"]
folds[grepl("early", label),type2:="introduced-Virginia"]
folds[grepl("mid", label),type2:="introduced-Virginia"]

folds[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
folds<-folds[order(type2)]

#code below doesn't work because the x values are different for each CV run
folds.sum<-folds[,.(timepoint.min=min(y, na.rm=T),
                    timepoint.max=max(y, na.rm=T)),
                 .(label, x, pop, loc, type, type2)]

#need to make windows instead
windows<-data.table(start=seq(from=-1, to=4.9, by =0.05),
                    end=seq(from=-0.9, to=5, by =0.05))
windows[,index:=c(1:nrow(windows))]


ribbons<-foreach(i=windows$index,.combine="rbind", .errorhandling="remove")%do%{
  window.data<-folds[x>=(10^windows[i,start])&x<(10^windows[i,end])]
  folds.sum<-window.data[,.(timepoint.min=min(y, na.rm=T),
                            timepoint.max=max(y, na.rm=T)),
                         .(label, pop, loc, type, type2)]
  folds.sum[,window.start:=(10^windows[i,start])]
  return(folds.sum)
}

ribbons[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
ribbons<-ribbons[order(type2)]

smcpp.plot<-ggplot()+
  geom_ribbon(data=ribbons[type2%in%c("native", "introduced-Florida", "introduced-Virginia")], aes(x=window.start, ymin=timepoint.min, ymax=timepoint.max, fill=as.factor(type2), group=label), alpha=0.15)+
  geom_line(data=final[type2%in%c("native", "introduced-Florida", "introduced-Virginia")], aes(x=x, y=y, color=type2, group=label), size=1.5)+  
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +annotation_logticks()+
  labs(x="years before present", y="population size", color=NULL, title="4 knots, cubic")+
  scale_color_manual(values = friendly_pal("ito_seven")[c(1,4,6)])+
  scale_fill_manual(values = friendly_pal("ito_seven")[c(1,4,6)])+
  theme(legend.position = c(0.3, 0.1)) + guides(fill="none")

smcpp.plot

ggplot(final)+geom_line(aes(x=log10(x), y=log10(y),group=label, color=label))
