library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(foreach)
library(doMC)
registerDoMC(20)
library(scales)
library(vcfR)


ld<-fread("/scratch/perickso/private/ind_seq/popgen/ld_decay_NorthAmerica.txt")
ld.sum<-ld[,.(avg.ld=mean(ld, na.rm=T), max.ld=max(ld, na.rm=T), N=.N), .(focal.snp, focal.chr, focal.pos)]


win.bp <- 1e5
step.bp <- 1e5

wins<-foreach(chr.i=c(1:5),
              .combine="rbind", 
              .errorhandling="remove")%dopar%{
                
                tmp <- ld[focal.chr==paste0("Scaffold_",chr.i)]
                data.table(CHR=chr.i,
                           start=seq(from=1, to=max(tmp$focal.pos)-(win.bp), by=step.bp),
                           end=seq(from=0, to=max(tmp$focal.pos)-(win.bp), by=step.bp) + win.bp)
              }

wins[,index:=c(1:nrow(wins))]
setkey(wins,index)


#test if each window has any snps that are 100, 200, 300, 400, or 500 kb away with rsq > 0.75
ld.win<-foreach(window=wins$index, .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  if(window%%1000==0){
                    print(window)
                  }
                  ld.tmp<-ld[dist%in%(c(1e5, 2e5, 3e5, 4e5, 5e5)) & focal.chr==paste0("Scaffold_",wins[J(window),CHR]) & focal.pos>=wins[J(window), start] & focal.pos<=wins[J(window),end]]
                  data.table(index=window,
                             mean.ld=mean((ld.tmp$ld)^2, na.rm=T),
                             med.ld=median((ld.tmp$ld)^2, na.rm=T),
                             max.ld=max((ld.tmp$ld)^2, na.rm=T),
                             q95=quantile((ld.tmp$ld)^2, .95, na.rm=T),
                             window.ld=max(abs(ld.tmp$ld))>0.75,
                             nsnps=sum(abs(ld.tmp$ld)>0.75),
                             n.snp.win=nrow(ld.tmp))
                }

ld.win<-merge(ld.win, wins, by="index")

#make a figure of LD 
ld.win[,inv:=rleid(window.ld)]
#ld.win[,inv:=rleid(nsnps>5)]
#ld.win[,inv2:=nsnps>5]
#ld.win[,inv.nsnps:=rleid(nsnps>5)]

ld.inv.sum<-ld.win[,.(inv.start=min(start), inv.end=max(end)), .(CHR, inv, window.ld)]
#ld.inv.sum.nsnps<-ld.win[,.(inv.start=min(start), inv.end=max(end)), .(CHR, inv2, inv.nsnps)]

ld.inv.sum[,length:=inv.end-inv.start]


#convert inversions to a bed file format
largeinv<-ld.inv.sum[window.ld==T&length>1000000]
bed<-data.table(chr=paste0("Scaffold_", largeinv$CHR),
                start=largeinv$inv.start-1,
                stop=largeinv$inv.end-1)

write.table(bed, "/scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed", quote=F, sep="\t", col.names = F, row.names = F)



#read in bed file for plotting
inv<-fread("/scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed", header=F)
setnames(inv, c("chr", "start", "stop"))
ld[, chr:=gsub("Scaffold_", "Chr. ", focal.chr)]
ld.win[,chr:=paste0("Chr. ", CHR)]

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

ld[dist==1e5, dist_lab:="100 kb"]
ld[dist==2e5, dist_lab:="200 kb"]
ld[dist==3e5, dist_lab:="300 kb"]
ld[dist==4e5, dist_lab:="400 kb"]
ld[dist==5e5, dist_lab:="500 kb"]

inv[,chr:=gsub("Scaffold_", "Chr. ", chr)]

#bring in data from smoove

vcf <- read.vcfR("/scratch/perickso/private/ind_seq/sv/zap_all_called_sv.smoove.square.vcf.gz")

#get data from INFO column with has info about the structural variantsinfo
sv.info<-as.data.table(INFO2df(vcf))

#extract fixed data into a data.table

v<-as.data.table(getFIX(vcf))

sv.data<-cbind(v, sv.info)
sv.data[,POS:=as.numeric(POS)]
sv.data[,SVLEN:=as.numeric(SVLEN)]

inv.smoove<-sv.data[as.numeric(SVLEN)>100000&as.numeric(SU)>10, .(CHROM, POS, END, SVLEN)]

inv.smoove[,chr:=gsub("Scaffold_", "Chr. ", CHROM)]

#assuming that any SVs large enough to influence PCs will be at least 500 kb

a<- ggplot(ld[dist%in%c(1e5, 2e5, 3e5, 4e5, 5e5)])+geom_point(aes(x=focal.snp, y=ld^2, color=focal.chr), size=0.5)+
  facet_grid(dist_lab~chr, scales="free")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  labs(x="SNP #", y="LD")+
  guides(color="none")+
  geom_hline(yintercept=0.75, linetype="dashed", color="grey50", size=0.5)+
  scale_x_continuous(label=scientific, expand = c(0, 0), breaks=c(1e6, 2e6, 3e6, 4e6, 5e6))

b<-ggplot(ld.win)+geom_rect(data=inv, aes(xmin=start, xmax=stop, ymin=-Inf, ymax=Inf),fill="grey90")+
  geom_point(aes(x=start, y=as.integer(window.ld)), size=0.25)+
  geom_rect(data=inv.smoove, aes(xmin=POS, xmax=END, ymin=0.45, ymax=0.55), fill=friendly_pal("ito_seven")[6])+
  facet_grid(.~chr, scales="free") +
  scale_y_continuous(limits=c(-.25, 1.25), breaks=c(0,1))+
  labs(y="high LD\nwindow", x="Position (bp)")+
  guides(color="none")+
  #scale_color_manual(values = friendly_pal("ito_seven"))+
  scale_x_continuous(label=scientific, expand = c(0, 0), breaks=c(1.5e7))+
  theme(strip.text.x = element_blank())

jpeg("/scratch/perickso/private/ind_seq/popgen/plots/LD_figure.jpg", height=10, width=8,units="in", res=300 )
plot_grid(a,b, nrow=2, rel_heights = c(0.8, 0.2), align="hv", axis="lr", labels=c("A", "B"))
dev.off()





