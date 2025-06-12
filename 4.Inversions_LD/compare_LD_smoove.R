library(SNPRelate)
library(vcfR)
library(data.table)

vcf <- read.vcfR("/scratch/perickso/private/ind_seq/sv/zap_all_called_sv.smoove.square.vcf.gz")

#get data from INFO column with has info about the structural variantsinfo
sv.info<-as.data.table(INFO2df(vcf))

#extract fixed data into a data.table

v<-as.data.table(getFIX(vcf))

sv.data<-cbind(v, sv.info)
sv.data[,POS:=as.numeric(POS)]
sv.data[,SVLEN:=as.numeric(SVLEN)]

#assuming that any SVs large enough to influence PCs will be at least 500 kb
inv<-sv.data[as.numeric(SVLEN)>1000000&as.numeric(SU)>10, .(CHROM, POS, END, SVLEN)]
invLD<-fread("/scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed", header=F)
sv.data[SVTYPE=="INV"&as.numeric(SVLEN)>1000000]
