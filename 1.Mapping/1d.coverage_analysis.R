library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(viridis)
library(ggsci)
library(ggpubfigs)

#list of samples
samps<-fread("/scratch/perickso/private/ind_seq/all_individual_samples.txt", header=F)$V1
samps2<-fread("/scratch/perickso/private/raw_data/published/SRR_Acc_List.txt" , header=F)$V1
#read in coverage per chromosome
cov<-foreach(i=c(samps,samps2), .errorhandling="remove")%do%{
    j<-fread(paste0("/scratch/perickso/private/ind_seq/coverage/", i, ".coverage.txt"))
    j[,sample.id:=i]
    return(j)
}
cov<-rbindlist(cov)
setnames(cov, "#rname", "scaffold")
cov[,chr:=as.numeric(tstrsplit(scaffold, split="_")[[2]])]

#pull out chromosome lengths
intervals<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)[,1:2]
names(intervals)=c("scaffold", "length")
intervals[,chr:=as.numeric(tstrsplit(scaffold, split="_")[[2]])]



#get coverage per chromosome across samples
cov.sum<-cov[,.(mean.cov=mean(meandepth, na.rm=T), med.cov=median(meandepth,na.rm=T), sd.cov=sd(meandepth,na.rm=T)), .(scaffold, endpos)]

hist(cov.sum$mean.cov)


#read in sample information
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated.csv", header=T)
setnames(metadata, "sample.name", "sample.id")

cov<-merge(cov, metadata, by="sample.id")

jpeg("/scratch/perickso/private/ind_seq/Figures/depth_by_scaffold.jpeg", height=6, width=8, units="in", res=300)
ggplot(cov[chr<=10&(Sex=="M"|Sex=="F")])+ 
  geom_boxplot(aes(x=as.factor(chr), y=meandepth, color=Sex))+
  labs(x="Scaffold", y="Mean Sequencing Depth")+
  geom_text(data=intervals[chr<=10], aes(x=as.factor(chr), label=as.character(length), y=42),angle = 45)+
  scale_color_manual(values = friendly_pal("bright_seven"))+
  lims(y=c(0,45))
  
dev.off()



