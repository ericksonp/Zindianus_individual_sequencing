library(data.table)
library(ggplot2)
library(ggpubfigs)
library(ggbio)
library(cowplot)
library(gggenomes)
theme_set(theme_cowplot())
library(foreach)

depth<-fread("/scratch/perickso/private/ind_seq/sv/depth_analysis_5kb.csv")
metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", header=T, drop=1)

depth<-merge(depth, metadata, by="sample.id")



win<-depth[loc.spec%in%c("Africa", "FL", "VA-CM", "VA-HPO")&chr==3&win>350000&win<1100000]
win[loc.spec=="VA-HPO"|loc.spec=="VA-CM", loc.spec:="VA"]
depth.sum<-win[,.(mean.depth=mean(rel_depth),
                    sd.depth=sd(rel_depth)),
                 .(loc.spec, chr, win)]


#a<-ggplot(win)+geom_line(aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.25)+
#  scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
#  facet_wrap(~assigned_sex)

#this part didn't work very well
# gff<-"/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff"
# # 
# ann<-read_feats(gff)
# ann.dt<-as.data.table(ann)
# # 
# b<-gggenomes(ann.dt[seq_id=="Scaffold_3"&end>600000&start<725000&type=="gene"])+geom_gene(size=10) +geom_seq_label() 
#  # scale_x_continuous(limits=c(600000,725000), expand = c(0, 0))
# 
# 
# plot_grid(a,b, nrow=2, align="v", axis="lr")
# # 
# candidates<-ann.dt[seq_id=="Scaffold_3"&end>600000&start<725000&type=="gene", feat_id]
# 
# genes<-fread("//scratch/perickso/private/ref/dovetail_transcripts_blast_protein_merged_beste.txt")
# genes[,gene:=tstrsplit(query, split="-")[[1]]]
# 
# genes[gene%in%candidates]

#what if we add in the repeat masker elements
repeats<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.RepeatMasked.gff.gz")
setnames(repeats, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

#ggplot()+geom_line(data=win, aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.25)+
#  scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
#geom_rect(data=repeats[chr=="Scaffold_3" &end>600000&start<725000], aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1))

#can we add the genes at the top?

ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

# ggplot()+geom_line(data=win, aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.25)+
#   scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
#   geom_rect(data=repeats[chr=="Scaffold_3" &end>550000&start<900000], aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1))+
#   geom_rect(data=ann[chr=="Scaffold_3" &end>550000&start<900000], aes(xmin=start, xmax=end, ymin=2.3, ymax=2.4), fill=friendly_pal("ito_seven")[3])

#can we add in a line for FST for reference?

pixy.fst<-foreach(i=c(1:5))%do%{
  if(i==3) {
    data<-fread("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_3_females_fst.txt")
  } else {
    data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_4pops_5kb_Scaffold_", i, "_fst.txt"))
  }
  return(data)
}
pixy.fst<-rbindlist(pixy.fst)
pixy.fst[,window.mid:=(window_pos_1+window_pos_2)/2]
pixy.fst[,index:=rleid(chromosome, window_pos_1)]
pixy.fst[,comparison:=paste(pop1, pop2, sep="-")]
pixy.fst


# ggplot()+geom_line(data=win, aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.25)+
#   scale_color_manual(values = friendly_pal("ito_seven")[4:6]) +
#   geom_rect(data=repeats[chr=="Scaffold_3" &end>550000&start<900000], 
#             aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1))+
#   geom_rect(data=ann[info=="gene"& chr=="Scaffold_3" &end>550000&start<900000], 
#             aes(xmin=start, xmax=end, ymin=2.9, ymax=3), fill=friendly_pal("ito_seven")[1])+
#   geom_line(data=pixy.fst[comparison=="FL-VA"&chromosome=="Scaffold_3"&window_pos_1>550000&window_pos_2<900000], 
#             aes(x=window.mid, y=avg_wc_fst), linewidth=2, color=friendly_pal("ito_seven")[3])+
#   labs(x="Chromosome 3 position", y="relative sequencing depth/FST", color="Population")+
#   scale_x_continuous(labels = ~ format(.x, scientific = FALSE))
#   

#can we add in the high IHS SNPs?

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
ihs.VA[,pop:="Virginia"]
ihs.VA[order(IHS)] #973443

ihs.high<-ihs.VA[IHS>5&CHR=="Scaffold_3"&POSITION<1500000, POSITION]
snps.to.plot<-data.table(x=ihs.high,
                         y=0.8)

 win<-win[order(loc.spec)]
# a<-ggplot()+
#   geom_rect(aes(xmin=600000, xmax=720000, ymin=-Inf, ymax=Inf), fill="grey90")+
#   geom_line(data=win[loc.spec=="VA"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[6])+
#   geom_line(data=win[loc.spec=="FL"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[5])+
#   geom_line(data=win[loc.spec=="Africa"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[4])+
#   geom_rect(data=repeats[chr=="Scaffold_3" &end>350000&start<1100000], 
#             aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1))+
#   geom_rect(data=ann[info=="gene"& chr=="Scaffold_3" &end>350000&start<1100000], 
#             aes(xmin=start, xmax=end, ymin=2.9, ymax=3), fill=friendly_pal("ito_seven")[1])+
#   geom_line(data=pixy.fst[comparison=="FL-VA"&chromosome=="Scaffold_3"&window_pos_1>350000&window_pos_2<1100000], 
#             aes(x=window.mid, y=avg_wc_fst), linewidth=2, color=friendly_pal("ito_seven")[3])+
#   labs(x="Chromosome 3 position", y="relative sequencing depth/FST", color="Population")+
#   scale_x_continuous(labels = ~ format(.x, scientific = FALSE))+


a<-ggplot()+
  geom_rect(aes(xmin=600000, xmax=720000, ymin=-Inf, ymax=Inf), fill="grey90")+
  geom_line(data=depth.sum[loc.spec=="VA"], aes(x=win, y=mean.depth, group=loc.spec, color=loc.spec, color="VA"))+
  geom_line(data=depth.sum[loc.spec=="FL"], aes(x=win, y=mean.depth, group=loc.spec, color=loc.spec, color="FL"))+
  geom_line(data=depth.sum[loc.spec=="Africa"], aes(x=win, y=mean.depth, group=loc.spec, color=loc.spec, color="Africa"))+
  geom_line(data=pixy.fst[comparison=="FL-VA"&chromosome=="Scaffold_3"&window_pos_1>350000&window_pos_2<1100000], 
            aes(x=window.mid, y=avg_wc_fst, color="FST"),linetype="dashed")+
  geom_point(data=snps.to.plot, aes(x=x, y=2.7, color="IHS>5"), alpha=0.25, stroke=NA)+
  scale_color_manual(values = c( friendly_pal("ito_seven")[3:6], "black"), breaks=c("IHS>5", "Africa", "FL", "VA", "FST")) +
  geom_rect(data=repeats[chr=="Scaffold_3" &end>350000&start<1100000], 
            aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1, fill="repeats"))+
  geom_rect(data=ann[info=="gene"& chr=="Scaffold_3" &end>350000&start<1100000], 
            aes(xmin=start, xmax=end, ymin=2.5, ymax=2.6, fill="genes"))+
  scale_fill_manual(values=c(friendly_pal("ito_seven")[1], "grey30"), breaks=c("genes", "repeats"))+
  labs(x=NULL, y="relative sequencing\ndepth/FST", color=NULL, fill=NULL)+
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE))



b<-ggplot()+
  geom_line(data=win[loc.spec=="VA"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[6])+
  geom_line(data=win[loc.spec=="FL"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[5])+
  geom_line(data=win[loc.spec=="Africa"], aes(x=win, y=rel_depth, group=sample.id, color=loc.spec), linewidth=0.15, color=friendly_pal("ito_seven")[4])+
  geom_point(data=snps.to.plot, aes(x=x, y=2.8), alpha=0.25, stroke=NA, color=friendly_pal("ito_seven")[3])+
  geom_rect(data=repeats[chr=="Scaffold_3" &end>350000&start<1100000], 
            aes(xmin=start, xmax=end, ymin=-0.2, ymax=-0.1))+
  geom_rect(data=ann[info=="gene"& chr=="Scaffold_3" &end>350000&start<1100000], 
            aes(xmin=start, xmax=end, ymin=2.6, ymax=2.7), fill=friendly_pal("ito_seven")[1])+
  geom_line(data=pixy.fst[comparison=="FL-VA"&chromosome=="Scaffold_3"&window_pos_1>350000&window_pos_2<1100000], 
            aes(x=window.mid, y=avg_wc_fst), linetype="dashed", color="black")+
  labs(x="Chromosome 3 position", y="relative sequencing\ndepth/FST", color="Population")+
  scale_x_continuous(labels = ~ format(.x, scientific = FALSE), limits=c(600000,720000))+
  theme(legend.position="none")


jpeg("/scratch/perickso/private/ind_seq/popgen/plots/scaffold_3_depth_individual.jpg", 
     height=5, 
     width=8,
     units="in", 
     res=300 )
plot_grid(a,b,nrow=2, labels=c("A", "B"), align="v", axis="lr")
dev.off()

#look at genotypes in this region--are africa all hets??
library(SNPRelate)
geno <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gds" , allow.fork=T)
a<-snpgdsSNPList(geno)
info<-data.table(snp.id=a$snp.id,
                 chr=tstrsplit(a$chromosome, split="_")[[2]],
                 pos=a$pos,
                 freq=a$afreq)
setkey(info, chr, pos)

snps.to.use<-info[chr==3&pos>600000&pos<725000, snp.id]

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)
metadata[loc.spec=="VA-CM"|loc.spec=="VA-HPO", loc.spec:="VA"]
samps.to.use<-metadata[assigned_sex=="F", sample.id]

genos<-snpgdsGetGeno(geno, snp.id=snps.to.use, sample.id=samps.to.use, with.id=T)

geno.dt<-as.data.table(genos$genotype)
names(geno.dt)=paste0("snp", snps.to.use)
geno.dt[,sample.id:=genos$sample.id]

genos.melt<-melt(geno.dt, id.var="sample.id", measure.var=names(geno.dt)[1:(ncol(geno.dt)-1)])
genos.melt[,snp.id:=as.integer(tstrsplit(variable, split="snp")[[2]])]
genos.melt<-merge(genos.melt, metadata, by="sample.id")
genos.melt<-genos.melt[loc.spec%in%c("Africa", "FL", "VA")][order(loc.spec)]
genos.melt[,ind.id:=rleid(sample.id)]
genos.melt<-genos.melt[order(snp.id)]
genos.melt[,pos.id:=rleid(snp.id)]
genos.melt<-merge(genos.melt, info, by="snp.id")


haplo.plot<-ggplot(genos.melt[loc.spec%in%c("Africa", "FL", "VA")])+
  geom_point(aes(x=pos, y=ind.id, color=as.factor(value)), shape=15)+
  facet_grid(loc.spec~., scales="free_y", space="free_y")+
  scale_color_manual(values=friendly_pal("ito_seven")[c(1,4,7)])
  

