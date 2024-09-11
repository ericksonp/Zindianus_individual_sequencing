library(data.table)
library(SNPRelate)
library(vcfR)

metadata=fread("/scratch/perickso/private/ind_seq/zap_full_info_updated.csv", header=T)
sex<-fread("/scratch/perickso/private/ind_seq/assigned_sex_by_coverage.txt")
setnames(sex, "sample.name", "sample.id")




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

metadata<-merge(metadata, sex, by="sample.id")
write.csv(metadata, "/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv")


#snpgdsVCF2GDS(vcf.fn = "/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gz",  
#              out.fn="/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gds", 
#              method="biallelic.only")


genofile <- snpgdsOpen("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gds" , allow.fork=T)
samps <- read.gdsn(index.gdsn(genofile, "sample.id")) 

samps.africa<-metadata[(sample.id%in%samps)&continent=="Africa"]$sample.id


a<-snpgdsSNPList(genofile)
info<-data.table(snp.id=a$snp.id,
                 chr=a$chromosome,
                 pos=a$pos)

b<-snpgdsSNPRateFreq(genofile, sample.id=samps.africa, with.snp.id=T)

#are these ref or alt freqs?
#snpgdsGetGeno(genofile,sample.id=samps.africa, snp.id=1)
#most genotypes are 2 (ref/ref), and allele frquency is 0.95, so allele freq is ref freq
#now get alt and reference alleles


africa.freqs<-data.table(snp.id=b$snp.id,
                         freq=b$AlleleFreq)

africa.freqs<-merge(africa.freqs, info, by="snp.id")

#now use bcftools to get metadata 

#bcftools query -f '%CHROM %POS  %REF  %ALT \n' /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gz > /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.snpinfo.txt


snp.data<-fread("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.snpinfo.txt", header=F)
names(snp.data)=c("chr", "pos", "REF", "ALT")

africa.freqs<-merge(africa.freqs, snp.data, by=c("chr", "pos"))

africa.freqs[,AA:=ifelse(freq>0.5, REF, ALT)]

write.table(africa.freqs[,.(chr, pos, AA)], file="/scratch/perickso/private/ind_seq/popgen/africa_ancestral_alleles.txt", quote=F, row.names=F, col.names=F, sep="\t")




#use this file to make something for plink to recode ancestral alleles
library(data.table)

africa.freqs<-fread("/scratch/perickso/private/ind_seq/popgen/africa_ancestral_alleles.txt")
setnames(africa.freqs, c("chr", "pos", "ancestral"))
bim<-fread("zaprionus.individual.2023.plink.bim", header=F )
setnames(bim, c("chr", "id", "blank", "pos", "ref", "alt"))

useme<-merge(africa.freqs, bim, by=c("chr", "pos"))

write.table(useme[,.(id, ancestral)], file="/scratch/perickso/private/ind_seq/popgen/new_reference_alleles_for_plink.txt", quote=F, sep="\t", row.names=F, col.names=F)
